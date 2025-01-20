import json
from fragment import FragmentList

import json

def read_json(file_path):
    """
    Reads a JSON file and processes it into fragment data and allow_overlapping flag.
    :param file_path: Path to the JSON file.
    :return: A tuple (fragments, allow_overlapping).
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
        
        # Extract fragments and allow_overlapping
        fragments = data.get("fragments", [])
        allow_overlapping = data.get("allow_overlapping", False)
        return fragments, allow_overlapping
    except FileNotFoundError:
        raise FileNotFoundError(f"File '{file_path}' not found.")
    except json.JSONDecodeError as e:
        raise ValueError(f"Failed to decode JSON: {e}")


def print_json(data):
    """
    Pretty-prints a Python dictionary as JSON.
    :param data: Dictionary to print.
    """
    import json
    print("JSON Data:")
    print(json.dumps(data, indent=4))


def print_fragment_list(fragment_list, allow_overlapping):
    """
    Prints the FragmentList and allow_overlapping status.
    :param fragment_list: FragmentList object.
    :param allow_overlapping: Boolean value indicating if overlapping is allowed.
    """
    print("Fragment List:")
    print(fragment_list)
    print(f"Allow Overlapping: {allow_overlapping}")

