import subprocess
import numpy as np

from io import StringIO


def argument_string(names : list,
                    values : list) -> str:
    """
    Function Details
    ================
    Create argument string for S4 variables.
    
    Variable names and parameters must be in corresponding order within arrays.
    
    Parameters
    ----------
    names, values: lists
        List of array variable names matching lua script names, list of matching
        values for argument names.
    
    Returns
    -------
    arg_string : string
        Argument string required for lua to use pcall.
    
    See Also
    --------
    S4_RCWA

    Notes
    -----
    Lua pcall requires the variable names and variable values to be in a string
    format. Ensure the names and values are ordered to not mis-value something.

    Example
    -------
    None

    ----------------------------------------------------------------------------
    Update History
    --------------

    25/05/2024
    ----------
    Documentation updated.

    """
    arguments = [
        f'{name} = {value}'
        for name, value in zip(names, values)]
    arg_string = ('; ').join(arguments)
    return arg_string


def S4_RCWA(lua_script : str,
            argument_string : str) -> str:
    """
    Function Details
    ================
    Run S4 RCWA simulation with desired parameters.
    
    Use the lua script and argument string set earlier to run S4.
    
    Parameters
    ----------
    lua_script, argument_string : string
        Lua script file name and arguments as a string.
    
    Returns
    -------
    process : string
        Wavelength, Transmission, Reflection output from lua script as a string.
    
    See Also
    --------
    subprocess
    argument_string
    read_S4_output

    Notes
    -----
    The lua script is usually set to print the output of the simulation to the
    terminal, so python grabs this output and outputs it as a string.

    Example
    -------
    None

    ----------------------------------------------------------------------------
    Update History
    ==============

    25/05/2024
    ----------
    Documentation updated.

    """
    command = f'S4 -a "{argument_string}" {lua_script}'
    process = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE)
    return process


def read_S4_output(process_string : str) -> list:
    """
    Function Details
    ================
    Read S4 string output.

    Read S4 output from subprocess stdout.

    Parameters
    ----------
    process_string : string
        Wavelength, Transmission, Reflection output from lua script as a string.
    
    Returns
    -------
    wavelength, transmission, reflection : list
        Wavelength, transmission, and reflection values.

    See Also
    --------
    S4_RCWA

    Notes
    -----
    None

    Example
    -------
    None

    ----------------------------------------------------------------------------
    Update History
    ==============

    25/05/2024
    ----------
    Documentation updated.

    """
    try:
        wavelength, transmission, reflection = np.genfromtxt(
            fname=StringIO(
                process_string.stdout.decode('utf-8')),
            delimiter='\t',
            unpack=True)
    except:
        wavelength, transmission, reflection = np.genfromtxt(
            fname=StringIO(
                process_string.stdout.decode('utf-8')),
            delimiter=',',
            unpack=True)
    return wavelength, transmission, reflection