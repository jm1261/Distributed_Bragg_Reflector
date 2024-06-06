import json
import numpy as np


def load_json(file_path: str) -> dict:
    """
    Loads .json file types.

    Use json python library to load a .json file.

    Parameters
    ----------
    file_path : string
        Path to file.

    Returns
    -------
    json file : dictionary
        .json dictionary file.

    See Also
    --------
    save_json_dicts

    Notes
    -----
    json files are typically dictionaries, as such the function is intended for
    use with dictionaries stored in .json file types.

    Example
    -------
    >>> my_dictionary = load_json(file_path="/Path/To/File")
    >>> my_dictionary
    {
        "Key 1": Value1,
        "Key 2": Value2
    }

    """
    with open(file_path, 'r') as file:
        return json.load(file)


def convert(o: str) -> TypeError:
    """
    Check data type.

    Check type of data string.

    Parameters
    ----------
    o : string
        String to check.

    Returns
    -------
    TypeError : Boolean
        TypeError if string is not suitable.

    See Also
    --------
    None.

    Notes
    -----
    None.

    Example
    -------
    None.

    """
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError


def save_json_dicts(out_path: str,
                    dictionary: dict) -> None:
    """
    Save .json file types.

    Use json python library to save a dictionary to a .json file.

    Parameters
    ----------
    out_path : string
        Path to file.
    dictionary : dictionary
        Dictionary to save.
    
    Returns
    -------
    None

    See Also
    --------
    load_json

    Notes
    -----
    json files are typically dictionaries, as such the function is intended for
    use with dictionaries stored in .json file types.

    Example
    -------
    >>> save_json_dicts(
        out_path="/Path/To/File",
        dictionary=my_dictionary)

    """
    with open(out_path, 'w') as outfile:
        json.dump(
            dictionary,
            outfile,
            indent=4,
            default=convert)
        outfile.write('\n')


def DBR_log_out(out_path : str,
                data : list):
    """
    Function Details
    ================
    Save log data out of the simulation.

    Parameters
    ----------
    out_path: str
        Path to save.
    data: list
        Data array to save, list of lists [[], [], []]

    Returns
    -------
    None.

    See Also
    --------
    None.

    Notes
    -----
    None.

    Example
    -------
    None.

    ----------------------------------------------------------------------------
    Update History
    ==============

    29/05/2024
    ----------
    Created.

    """
    header_string = (
        f'Wavelength [nm]\tk0\tkc\tks\t'
        f'Transmission Coefficient\tReflection Coefficient\tPhase Coefficient')
    formatted_data = np.array(data)
    np.savetxt(
        out_path,
        formatted_data,
        header=header_string,
        delimiter='\t',
        fmt='%e')
