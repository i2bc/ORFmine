"""
Script utility used for the user to update the external softwares (IUPRED and/or Tango) local path
"""
import configparser
from pathlib import Path
import sys

from orfmine.orfold.lib import utils
from orfmine import ROOT_PATH

CONFIG_FILE = ROOT_PATH / 'softwares.ini'


def get_binding_choice():
    key_input = ""
    binding_type_choice = ""
    is_valid_number = False

    print("You are going to inform ORFmine where external softwares such as IUPred2a or Tango are locally installed on your computer.")
    print("You can leave this prompt at any moment by entering either 'q', 'quit' or 'exit'.\n")

    while not isinstance(binding_type_choice, int) and not is_valid_number:
        print("Do you want to bind IUPred2a (1), Tango (2), or both (3)?:")
        key_input = input(" - Please choose the number that best suits your needs: ")
        if key_input in ['q', 'quit' or 'exit']:
            break
        try:
            binding_type_choice = int(key_input)
            if not 1 <= binding_type_choice <= 3:
                print("\nSorry, the number must be either 1, 2 or 3, you entered {}...\n".format(key_input))
                binding_type_choice = key_input
            else:
                is_valid_number = True            
        except:
            print("\nSorry, you must enter a number (either 1, 2 or 3)...\n")

    return binding_type_choice


def check_path(path: str="", software: str="iupred"):
    software_name = "IUPred2a" if software == "iupred" else "Tango"
    error_message = ""
    missing_files = []
    is_valid = True

    exec_sources = {
        "iupred": {
            "darwin": ["iupred2a_lib.py", "iupred2a.py", "data"],
            "win32": ["iupred2a_lib.py", "iupred2a.py", "data"],
            "linux": ["iupred2a_lib.py", "iupred2a.py", "data"],
        },
        "tango": {
            "darwin": ["tango2_3_1"],
            "win32": ["Tango.exe"],
            "linux": ["tango_x86_64_release"],
        }
    }

    # case where the path given by the user do not exist   
    if not Path(path).exists():
        error_message = "\nSorry, it looks like {} do not exist... Please ensure to give the correct absolute path where {} resides.".format(path, software_name)
        is_valid = False

        return is_valid, error_message

    # get expected source files according to the software and platform
    if sys.platform not in exec_sources["iupred"]:
        source_files = exec_sources[software]["linux"]
    else:
        source_files = exec_sources[software][sys.platform]

    # cases where the given path exists but expected at least one expected source files is missing 
    for filename in source_files:
        source_file = Path(path) / filename
        if not source_file.exists():
            missing_files.append(str(source_file))
            is_valid = False

    if missing_files:        
        term = "files" if len(missing_files) > 1 else "file"        
        error_message = "\nSorry, it looks like the following {} do not exist in {}: {}".format(term, path, ", ".join(missing_files))    

    return is_valid, error_message


def prompt_user():
    key_input = ""
    iupred_path = ""
    tango_path = ""
    is_iupred_path_valid = False
    is_tango_path_valid = False

    binding_type_choice = get_binding_choice()

    if binding_type_choice == 1:
        while not is_iupred_path_valid:
            print("\nPlease enter the absolute path where IUPred2a is installed on your computer.")
            key_input = input(" - Absolute path to IUPred2A (e.g. /home/user/iupred2a/): ")
            if key_input in ["q", "quit", "exit"]:
                break

            is_iupred_path_valid, error_message = check_path(path=key_input, software='iupred')
            if error_message:
                print(error_message)
            else:
                iupred_path = key_input


    elif binding_type_choice == 2:
        while not is_tango_path_valid:
            print("\nPlease enter the absolute path where Tango is installed on your computer.")
            key_input = input(" - Absolute path to Tango (e.g. /home/user/tango/): ")
            if key_input in ["q", "quit", "exit"]:
                break

            is_tango_path_valid, error_message = check_path(path=key_input, software='tango')
            if error_message:
                print(error_message)                
            else:
                tango_path = key_input

    elif binding_type_choice == 3:
        while not (is_iupred_path_valid and is_tango_path_valid):
            if not is_iupred_path_valid:
                while not is_iupred_path_valid:
                    print("\nPlease enter the absolute path where IUPred2a is installed on your computer.")
                    key_input = input(" - Absolute path to IUPred2A (e.g. /home/user/iupred2a/): ")
                    if key_input in ["q", "quit", "exit"]:
                        break

                    is_iupred_path_valid, error_message = check_path(path=key_input, software='iupred')
                    if error_message:
                        print(error_message)
                    else:
                        iupred_path = key_input

            if not is_tango_path_valid:                
                while not is_tango_path_valid:
                    print("\nPlease enter the absolute path where Tango is installed on your computer.")
                    key_input = input(" - Absolute path to Tango (e.g. /home/user/tango/): ")
                    if key_input in ["q", "quit", "exit"]:
                        break

                    is_tango_path_valid, error_message = check_path(path=key_input, software='tango')
                    if error_message:
                        print(error_message)                
                    else:
                        tango_path = key_input

    return {"iupred": iupred_path, "tango": tango_path}

"""
def update_config_file():
    response = prompt_user()

    config_comments = [
        '; you can enter the absolute root path where IUPred2a and/or Tango softwares executables resides',
        '; examples',
        '; iupred = "/home/user/iupred2a/"',
        '; tango = "/home/user/tango/"'
    ]

    # retrieve current data in config.ini
    #current_external_softwares = utils.read_config_file()

    # start a new config instance 
    config = configparser.ConfigParser(allow_no_value=True)
    
    config['EXTERNAL_SOFTWARE'] = {}
    for comment in config_comments:
        config['EXTERNAL_SOFTWARE'][comment] = None

    for key, value in response.items():
        if value:
            config['EXTERNAL_SOFTWARE'][key] = '"{}"'.format(value)
        else:
            config['EXTERNAL_SOFTWARE'][key] = current_external_softwares[key]
    
    with open(CONFIG_FILE, 'w') as configfile:
        config.write(configfile)
"""
