import subprocess
import numpy as np

from io import StringIO
from src.dataIO import DBR_log_out


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


def create_Bragg_stacks(first_index : float,
                        first_thickness : float,
                        second_index : float,
                        second_thickness : float,
                        cavity_location : int,
                        cavity_index : float,
                        cavity_thickness : float,
                        number_of_pairs : int,
                        symmetry : str,
                        extra_layer : str) -> tuple:
    """
    Function Details
    ================
    Build a Bragg mirror stack with the input refractive indices.

    Parameters
    ----------
    first_index, first_thickness: float
        Refractive index (RIU) and layer thickness (in nm) for the first layer.
    second_index, second_thickness: float
        Refractive index (RIU) and layer thickness (in nm) for the second layer.
    cavity_index, cavity_thickness: float
        Refractive index (RIU) and layer thickness (in nm) for the cavity layer.
    cavity_location, number_of_pairs: integer
        Cavity location in the Bragg stack as an array index, remember Python
        starts at index 0. Number of mirror pairs, note that the cavity replaces
        one layer of the mirror pair at the given index, eg., if there are 4
        pairs, and you put cavity at index=2, it will go (first, second, cavity,
        second).
    symmetry, extra_layer: string
        If "True" will ensure the mirror pairs are symmetrical about the cavity
        layer, else will insert first, second, cavity, first, second, etc. If
        "True" will add an extra mirror pair layer (singular) underneath the
        final mirror pair to ensure symmetry.

    Returns
    -------
    n, t: list
        Refractive index array, thickness array.
    number_of_layers : int
        Number of layers, given as the length of n.

    See Also
    --------
    None.

    Notes
    -----
    Note that allowing for symmetry will include an extra layer on the bottom,
    i.e., it will allow for the mirrors to be placed either side of the stack.

    Example
    -------
    None.

    ----------------------------------------------------------------------------
    Update History
    ==============

    29/05/2024
    ----------
    Created from previous script.

    """
    n = [first_index, second_index] * number_of_pairs
    t = [first_thickness, second_thickness] * number_of_pairs
    if symmetry == "True":
        n[cavity_location] = cavity_index
        t[cavity_location] = cavity_thickness
    else:
        n[cavity_location] = cavity_index
        t[cavity_location] = cavity_thickness / 2
        n[cavity_location + 1] = cavity_index
        t[cavity_location + 1] = cavity_thickness / 2
    if extra_layer == "True":
        n.append(first_index)
        t.append(first_thickness)
    number_of_layers = len(n)
    return n, t, number_of_layers


def Bragg_transfer_matrices(wavelength_range : list,
                            number_of_layers : int,
                            stack_n : list,
                            stack_t : list,
                            effective_index : float,
                            cladding_index : float,
                            substrate_index : float,
                            out_path : str):
    """
    Function Details
    ================
    Perform transfer matrix calculation for the Bragg Stack.

    Parameters
    ----------
    wavelength_range, stack_n, stack_t: list
        Wavelength range at which to analyse the stack. Refractive index list
        for the Bragg stack. Thickness list for the Bragg stack.
    number_of_layers: int
        Number of layers in the stack, note that this is double the number of
        mirror pairs.
    effective_index, cladding_index, substrate_index: float
        Refractive index (RIU) values for the effective mode index, top layer
        index, and the substrate layer index.
    out_path: string
        Path to save log file.

    Returns
    -------
    transmission, reflection, phase: list
        Transmission, reflection, and phase results through the Bragg stack for
        a given wavelength range and parameters.

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
    Created from script.

    """
    transmission = []
    reflection = []
    phase = []

    simulation_log = []

    for wavelength in wavelength_range:
        k_0 = (2 * np.pi) / wavelength
        M_0 = [[1, 0], [0, 1]]

        for index in range(number_of_layers):
            refractive_index = ((stack_n[index] ** 2) - (effective_index ** 2))
            nk = np.sqrt(refractive_index) * k_0
            nkt = nk * stack_t[index]

            M = [
                [np.cos(nkt), (1j / nk) * np.sin(nkt)],
                [1j * nk * np.sin(nkt), np.cos(nkt)]
            ]
            M_0 = np.matmul(M_0, M)

        k_c = k_0 * np.sqrt((cladding_index ** 2) - (effective_index ** 2))
        k_s = k_0 * np.sqrt((substrate_index ** 2) - (effective_index ** 2))

        T = (
            (2 * k_s) /
            (
                (k_s * M_0[0][0]) +
                (k_c * M_0[1][1]) +
                (k_s * k_c * M_0[0][1]) +
                (M_0[1][0])
            )
        )
        transmit = T * T.conjugate()
        transmission.append(transmit)

        R = (
            (
                (k_s * M_0[0][0]) -
                (k_c * M_0[1][1]) +
                (k_s * k_c * M_0[0][1]) -
                (M_0[1][0])
            ) *
            (T / (2 * k_s))
        )
        reflect = R * R.conjugate()
        reflection.append(reflect)

        P = np.arctan(T.real / T.imag)
        phase.append(P)

        log_data = [
            wavelength,
            k_0,
            k_c,
            k_s,
            transmit.real,
            reflect.real,
            P]
        simulation_log.append(log_data)

    DBR_log_out(
        out_path=out_path,
        data=simulation_log)

    return transmission, reflection, phase
