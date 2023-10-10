from argparse import Namespace
import importlib.util
from pathlib import Path
import subprocess
import sys
from typing import Dict, List, Optional, Union


def get_package_path(package_name):
    """Retrieve the directory path of a package based on its name."""
    spec = importlib.util.find_spec(package_name)
    if spec and spec.origin:
        return str(Path(spec.origin).parent.resolve())
    
    return None


def add_container_args(parser):
    parser.add_argument(
        "--dry-run", "-D",
        required=False,
        action='store_true',
        default=False,
        help="Flag used to show the docker command line. Must be used in conjunction with '--docker' or '--singularity'"
    )

    parser.add_argument(
        "--dev",
        required=False,
        action='store_true',
        default=False,
        help=(
            "Enables development mode. Binds the local package directory to its "
            "corresponding path inside the container, allowing for real-time testing "
            "and changes without rebuilding the container. Should be used with either "
            "'--docker' or '--singularity'."
        )
    )

    container_group = parser.add_mutually_exclusive_group()
    container_group.add_argument("--docker", action='store_true', default=False, help="Flag used to run computations on a docker container")
    container_group.add_argument("--singularity", action='store_true', default=False, help="Flag used to run computations on a singularity container")

    return parser


class ContainerCLI:
    """A class to translate and run a command line interface (CLI) for a containerized application, either using Docker or Singularity.

    Attributes:
        args (argparse.Namespace): The arguments namespace from the command line.
        input_args (List[str]): A list of command line arguments related to input files.
        output_arg (str): The command line argument related to the output directory.
        image_base (str): Docker image fullname (e.g. f"{DOCKER_REPOSITORY}/{IMAGE_NAME}:{IMAGE_TAG}").
        cmd_args (Optional[List[str]], optional): Command line arguments. If not provided, sys.argv is used.
        prog (str, optional): Name of the program to be executed.
        software_bindings (dict, optional): Dictionary to hold software name and its path on the host system.
        container_type (str, optional): Type of container to be used ; either 'docker' or 'singularity'. Defaults to 'docker'.
        container_ipath (str, optional): The input path inside the container. Defaults to '/input'.
        container_opath (str, optional): The output path inside the container. Defaults to '/output'.
        extra_binding (dict, optional): A dictionary that maps external files to their respective paths inside the container for mounting.
        dev_mode (bool, optional): Flag used to set up a development mode.
        package_binding (dict, optional): A dictionary that maps a local package to its respective path inside the container.
    """

    VALID_CONTAINERS = ["docker", "singularity"]

    def __init__(
            self,
            args: Namespace,
            input_args: List[str],
            output_arg: str,
            image_base: str,
            cmd_args: Optional[List[str]] = None,
            prog: str="",
            workdir: str="",
            software_bindings: Optional[Dict[str, str]] = None,
            container_type: str="docker",
            container_ipath: str='/input',
            container_opath: str='/output',
            extra_bindings: dict={},
            package_binding: dict={},
            dev_mode: bool=False
                ) -> None:
        """

        Args:
            args (argparse.Namespace): The arguments namespace from the command line.
            input_args (List[str]): A list of command line arguments related to input files.
            output_arg (str): The command line argument related to the output directory.
            image_base (str): Docker image fullname (e.g. f"{DOCKER_REPOSITORY}/{IMAGE_NAME}:{IMAGE_TAG}").
            cmd_args (Optional[List[str]], optional): Command line arguments. If not provided, sys.argv is used.
            prog (str, optional): Name of the program to be executed.
            software_bindings (dict, optional): Dictionary to hold software name and its path on the host system.
            container_type (str, optional): Type of container to be used ; either 'docker' or 'singularity'. Defaults to 'docker'.
            container_ipath (str, optional): The input path inside the container. Defaults to '/input'.
            container_opath (str, optional): The output path inside the container. Defaults to '/output'.
            extra_binding (dict, optional): A dictionary that maps external files to their respective paths inside the container for mounting.
            dev_mode (bool, optional): Flag used to set up a development mode.
            package_binding (dict, optional): A dictionary that maps a local package to its respective path inside the container.
        """        
        self.args = args
        self.input_args = input_args
        self.output_arg = output_arg
        self.image_base = image_base

        self.cmd_args = cmd_args or sys.argv
        self.prog = prog or Path(cmd_args[0]).name
        self.workdir = str(Path(workdir).resolve()) if workdir else ""
        self.software_bindings = software_bindings
        
        self.container_type = container_type
        self.container_ipath = container_ipath
        self.container_opath = container_opath
        self.extra_bindings = extra_bindings

        self.package_binding = package_binding
        self.dev_mode = dev_mode

        self.container_handler = {
            'docker': {
                'base_cmd': ['docker', 'run', '-it', '--rm'],
                'binding_flag': '-v',
                'workdir_flag': '-w',
                'image_url': [self.image_base]
            },
            'singularity': {
                'base_cmd': ['singularity', 'exec'],
                'binding_flag': '-B',
                'workdir_flag': '--pwd',
                'image_url': [f"docker://{self.image_base}"]
            }
        }

        self._validate_container_type()
        self._get_io()
        self._create_outdir()
        self._generate_cli()

        
    def _validate_file_exists(self, filepath: Union[str, List[str]]) -> None:
        do_not_exist = []

        if not isinstance(filepath, list):
            filepath = [filepath]

        for _file in filepath:
            if not Path(_file).exists():
                do_not_exist.append(_file)

        if do_not_exist:
            raise FileNotFoundError(f"Error, the provided file(s) do(es) not exist: '{do_not_exist}'")

    def _validate_container_type(self) -> None:
        if self.container_type not in self.VALID_CONTAINERS:
            raise ValueError(f"Unsupported container type '{self.container_type}'. Supported types are: {', '.join(self.VALID_CONTAINERS)}")

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

        >>> self.cmd_args = ["prog", "-foo1", "/work/file1", "-foo2", "/work/file2", "-out", "/work/output"]
        >>> self.input_args = ["-foo1", "-foo2"]
        >>> self.output_arg = "-out"
        >>> self.iargs = ["foo1", "foo2"]
        >>> self.oarg = ["out"]
        >>> self.input_files = [["/work/file1"], ["/work/file2"]]
        >>> self.output = "/work/output"
        """
        iargs = [ self._arg_to_attr(arg=arg) for arg in self.input_args ]
        oarg = self._arg_to_attr(arg=self.output_arg)

        # retrieve attributes | ensure input files elements are typed as list for a more generic usage
        self.input_files = [[getattr(self.args, x)] if not isinstance(getattr(self.args, x), list) else getattr(self.args, x) for x in iargs]
        self.output = getattr(self.args, oarg)

        # Check if input files exist
        for file in self.input_files:
            self._validate_file_exists(filepath=file)

    def _create_outdir(self):
        Path(self.output).mkdir(exist_ok=True, parents=True)
            
    def _generate_cli(self):
        """Generate the expected command line from the given container type.
        """
        base_command = self.container_handler[self.container_type]['base_cmd']
        image_url = self.container_handler[self.container_type]['image_url']

        # set optional workdir path on the container
        workdir_flag = self.container_handler[self.container_type]['workdir_flag']
        if not self.dev_mode:
            wd = [workdir_flag, self.workdir] if self.workdir else []
        else:
            wd = [workdir_flag, self.container_ipath]

        # translate the command line
        self.cli = base_command + self._generate_volume_bindings() + wd + image_url

        if not self.dev_mode:
            self.cli +=  [self.prog] + self._parse_arguments()

    def get_inp_binding_path(self, file_path):
        return f"{Path(file_path).resolve()}:{self.container_ipath}/{Path(file_path).name}"
        
    def _generate_volume_bindings(self):
        """Generate the syntax to mount/bind required io to the container filesystem

        Returns:
            mnt_point (List): List of strings defining the mounting
        """
        binding_flag = self.container_handler[self.container_type]['binding_flag']
        
        mnt_point = []

        # add input paths bindings
        input_binding_path = [self.get_inp_binding_path(_file) for input_file in self.input_files for _file in input_file]
        mnt_point += [item for input_path in input_binding_path for item in [binding_flag, input_path]]
        
        # add output path binding
        mnt_point += [binding_flag, f"{Path(self.output).resolve()}:{self.container_opath}"]

        # add external extra bindings
        if self.extra_bindings:
            for host_file, container_path in self.extra_bindings.items():
                container_path = f"{container_path}/{Path(host_file).name}"
                mnt_point += [binding_flag, f"{Path(host_file).resolve()}:{container_path}"]

        # add package binding - dev-mode
        if self.dev_mode:
            if not self.package_binding:
                print("Error: package_binding must be given as arguments with dev-mode")
                exit(0)
            for package_name, container_path in self.package_binding.items():
                mnt_point += [binding_flag, f"{get_package_path(package_name=package_name)}:{container_path}"]

        return mnt_point

    def _parse_arguments(self):
        """Get all the arguments from the original command line and rewrite it with the expected container filesystem mounting 

        Returns:
            cmd_arguments (List): List of strings defining the container suited arguments
        """
        # add all command line arguments except --docker or --singularity if present
        cmd_arguments = [ argv for argv in self.cmd_args[1:] if argv not in [f"--{x}" for x in self.VALID_CONTAINERS+["dev-mode"] ] ]

        # flatten input_files
        input_files = [item for input_file in self.input_files for item in input_file]

        # edit arguments relative to input(s)
        for arg in cmd_arguments:
            if arg in input_files:
                cmd_arguments[cmd_arguments.index(arg)] = f"{self.container_ipath}/{Path(arg).name}"
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

        error_detail = None
        try:
            subprocess.run(self.cli, stdout=sys.stdout, text=True, check=True)
        except subprocess.CalledProcessError as e:
            error_detail = e.stderr or e.stdout or "No additional error message provided."
        except Exception as e:
            error_detail = str(e)

        if error_detail:
            final_error_msg = f"Command failed with error:\n{error_detail}"
            print(final_error_msg)
            exit(1)


    

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



