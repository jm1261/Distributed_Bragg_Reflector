import numpy as np
import src.dataIO as io
import src.filepaths as fp
import src.analysis as anal
import src.plotting as plot

from pathlib import Path


def distributed_bragg_reflector(simulation_parameters : dict,
                                out_path : str) -> None:
    """
    Function Details
    ================
    ----------------------------------------------------------------------------
    Update History
    ==============

    29/05/2024
    ----------
    Created.

    """
    file_stem = "DBR_Cavity_Python"
    graph_file = Path(f'{out_path}/{file_stem}_Test4.png')

    n, t, layers = anal.create_Bragg_stacks(
        first_index=simulation_parameters["first_n"],
        first_thickness=simulation_parameters["first_t"],
        second_index=simulation_parameters["second_n"],
        second_thickness=simulation_parameters["second_t"],
        cavity_location=simulation_parameters["cavity_index"],
        cavity_index=simulation_parameters["cavity_n"],
        cavity_thickness=simulation_parameters["cavity_t"],
        number_of_pairs=simulation_parameters["number_of_pairs"],
        symmetry=simulation_parameters["symmetrical"],
        extra_layer=simulation_parameters["extra_layer"])

    plot.stack_plots(
        stack_n=n,
        stack_t=t,
        plot_dict=plot_dict,
        out_path=Path(f'{out_path}/{file_stem}_IndexTest.png'))

    wavelength_range = np.arange(
        simulation_parameters["wavelength_i"],
        (
            simulation_parameters["wavelength_f"] +
            simulation_parameters["wavelength_s"]),
        simulation_parameters["wavelength_s"])

    transmission, reflection, phase = anal.Bragg_transfer_matrices(
        wavelength_range=wavelength_range,
        number_of_layers=layers,
        stack_n=n,
        stack_t=t,
        effective_index=simulation_parameters["n_effective"],
        cladding_index=simulation_parameters["cladding_n"],
        substrate_index=simulation_parameters["substrate_n"],
        out_path=Path(f'{out_path}/{file_stem}_Log.txt'))

    plot.Bragg_plots(
        wavelength=wavelength_range,
        transmission=transmission,
        reflection=reflection,
        phase=phase,
        out_path=graph_file,
        plot_dict=plot_dict)


if __name__ == '__main__':

    """ Organisation """
    root = Path(f'{Path().absolute()}/Distributed_Bragg_Reflector')
    results_path = Path(f'{root}/Results')
    plot_dict = io.load_json(
        file_path=Path(f'{root}/Plot_Dict.json'))
    sim_params = io.load_json(
        file_path=Path(f'{root}/Simulation_Parameters.json'))

    """ Outputs """
    out_path = Path(
        f'{results_path}/'
        f'nc_{sim_params["cavity_n"]}_'
        f'n1_{sim_params["first_n"]}_'
        f'n2_{sim_params["second_n"]}')
    fp.check_dir_exists(dir_path=out_path)

    distributed_bragg_reflector(
        simulation_parameters=sim_params,
        out_path=out_path)
