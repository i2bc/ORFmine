from argparse import Namespace
from pathlib import Path
import subprocess
import sys
from typing import List, Optional


class ContainerCLI:
    """A class to translate and run a command line interface (CLI) for a containerized application, either using Docker or Singularity.

    Attributes:
        args (argparse.Namespace): The arguments namespace from the command line.
        input_args (List[str]): A list of command line arguments related to input files.
        output_arg (str): The command line argument related to the output directory.
        image_base (str): Docker image fullname (e.g. f"{DOCKER_REPOSITORY}/{IMAGE_NAME}:{IMAGE_TAG}").
        cmd_args (Optional[List[str]], optional): Command line arguments. If not provided, sys.argv is used.
        prog (str, optional): Name of the program to be executed
        container_type (str, optional): Type of container to be used ; either 'docker' or 'singularity'. Defaults to 'docker'.
        container_ipath (str, optional): The input path inside the container. Defaults to '/input'.
        container_opath (str, optional): The output path inside the container. Defaults to '/output'.
    """

    VALID_CONTAINERS = ("--docker", "--singularity")

    def __init__(self, args: Namespace, input_args: List[str], output_arg: str, image_base: str, cmd_args: Optional[List[str]] = None, prog: str="", container_type: str="docker", container_ipath: str='/input', container_opath: str='/output') -> None:
        """

        Args:
            args (argparse.Namespace): The arguments namespace from the command line.
            input_args (List[str]): A list of command line arguments related to input files.
            output_arg (str): The command line argument related to the output directory.
            image_base (str): Docker image fullname (e.g. f"{DOCKER_REPOSITORY}/{IMAGE_NAME}:{IMAGE_TAG}").
            cmd_args (Optional[List[str]], optional): Command line arguments. If not provided, sys.argv is used.
            prog (str, optional): Name of the program to be executed.
            container_type (str, optional): Type of container to be used ; either 'docker' or 'singularity'. Defaults to 'docker'.
            container_ipath (str, optional): The input path inside the container. Defaults to '/input'.
            container_opath (str, optional): The output path inside the container. Defaults to '/output'.
        """        
        self.args = args
        self.input_args = input_args
        self.output_arg = output_arg
        self.image_base = image_base

        self.cmd_args = cmd_args or sys.argv
        self.prog = prog or Path(cmd_args[0]).name
        
        self.container_type = container_type
        self.container_ipath = container_ipath
        self.container_opath = container_opath

        self.container_handler = {
            'docker': {
                'base_cmd': ['docker', 'run', '-it', '--rm'],
                'binding_flag': '-v',
                'image_url': [self.image_base]
            },
            'singularity': {
                'base_cmd': ['singularity', 'exec'],
                'binding_flag': '-B',
                'image_url': [f"docker://{self.image_base}"]
            }
        }

        self._validate_container_type()
        self._get_io()
        self._generate_cli()

        
    def _validate_file_exists(self, filepath: str) -> None:
        if not Path(filepath).exists():
            raise FileNotFoundError(f"Provided file path '{filepath}' does not exist.")

    def _validate_container_type(self, container_type: str) -> None:
        if container_type not in self.VALID_CONTAINERS:
            raise ValueError(f"Unsupported container type '{container_type}'. Supported types are: {', '.join(self.VALID_CONTAINERS)}")

    def _arg_to_attr(self, arg: str) -> str:
        """Translates an argument name to an attribute name compatible with argparse.Namespace format

        Args:
            arg (str): Argument names

        Returns:
            attr (str): Attribute name

        Usage examples:
        >>> arg_to_attr("--flag-1")
        'flag_1'
        >>> arg_to_attr("-flag_2")
        'flag_2'
        """
        while arg.startswith("-"):
            arg = arg.replace("-", "", 1)
        attr = arg.replace("-", "_")

        return attr

    def _get_io(self):
        """Get inputs/output paths from their argument flags in Namespace
        """
        iargs = [ self._arg_to_attr(arg=arg) for arg in self.input_args ]
        oarg = self._arg_to_attr(arg=self.output_arg)

        self.input_files = [ getattr(self.args, x) for x in iargs ]
        self.output = getattr(self.args, oarg)

        # Check if input files exist
        for file in self.input_files:
            self._validate_file_exists(filepath=file)
            
    def _generate_cli(self):
        """Generate the expected command line from the given container type.
        """
        base_command = self.container_handler[self.container_type]['base_cmd']
        image_url = self.container_handler[self.container_type]['image_url']
        self.cli = base_command + self._generate_volume_bindings() + image_url + [self.prog] + self._parse_arguments()

    def _generate_volume_bindings(self):
        """Generate the syntax to mount/bind required io to the container filesystem

        Returns:
            mnt_point (List): List of strings defining the mounting
        """
        binding_flag = self.container_handler[self.container_type]['binding_flag']
        
        mnt_point = []
        for input_file in self.input_files:
            mnt_point += [binding_flag, f"{Path(input_file).resolve()}:{self.container_ipath}/{Path(input_file).name}"]
        
        mnt_point += [binding_flag, f"{Path(self.output).resolve()}:{self.container_opath}"]

        return mnt_point

    def _parse_arguments(self):
        """Get all the arguments from the original command line and rewrite it with the expected container filesystem mounting 

        Returns:
            cmd_arguments (List): List of strings defining the container suited arguments
        """
        # add all command line arguments except --docker or --singularity if present
        cmd_arguments = [ argv for argv in self.cmd_args[1:] if argv not in self.VALID_CONTAINERS ]

        # edit arguments relative to input(s)
        for arg in cmd_arguments:
            if arg in self.input_files:
                cmd_arguments[cmd_arguments.index(arg)] = f"/input/{Path(arg).name}"
            elif arg == self.output:
                cmd_arguments[cmd_arguments.index(arg)] = self.container_opath

        if self.output_arg not in cmd_arguments:
            cmd_arguments += [self.output_arg, self.container_opath]

        return cmd_arguments

    def cmd(self):
        """Return the converted command line argument in a string format

        Returns:
            (str): Converted command line argument in a string format
        """
        return " ".join(self.cli)

    def show(self):
        """Print the converted command line argument
        """
        print(f"\n{self.cmd()}\n")

    def run(self):
        """Run the converted command line argument"""
        Path(self.output).mkdir(exist_ok=True, parents=True)
        try:
            subprocess.run(self.cli, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Command '{self.cmd()}' failed with error:\n{e.stderr}")
        except Exception as e:
            raise RuntimeError(f"Error executing the command: {self.cmd()}. Error: {str(e)}")


if __name__ == "__main__":
    # orftrack cmd line: orftrack -fna data/foo.fna -gff data/foo.gff -out test
    # containerized cli: 
    # docker run --rm -it \
    # -v /home/nche/data/foo.fna:/input/foo.fna \
    # -v /home/nche/data/foo.gff:/input/foo.gff \
    # -v /home/nche/test:/output orfmine:latest \
    # orftrack -fna /input/foo.fna -gff /input/foo.gff -out /output
    cmd_string = "orftrack -fna data/foo.fna -gff data/foo.gff -out test --docker"
    sys.argv = cmd_string.split()

    args = Namespace(fna='data/foo.fna', gff='data/foo.gff', out='test')
    input_args = ["fna", "gff"]
    output_arg = "out"
    out_dirname = "test"
    cli = ContainerCLI(args=args, input_args=input_args, output_arg=output_arg)



