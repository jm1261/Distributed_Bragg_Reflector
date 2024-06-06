import numpy as np
import src.dataIO as io
import src.filepaths as fp
import src.analysis as anal
import src.plotting as plot

from pathlib import Path


def distributed_Bragg_reflector(simulation_parameters : dict,
                                out_path : str,
                                file_stem : str,
                                plot_dict : dict,
                                cavity=False) -> None:
    """
    Function Details
    ================
    Run a distributed Bragg reflector simulation using the transfer matrix
    method.

    Parameters
    ----------
    simulation_parameters, plot_dict: dictionary
        Simulation settings dictionary containing:
            {
                "loop": "loop settings",\n
                "n_effective": mirror effective index (loss),\n
                "cladding_n": cladding refractive index,\n
                "cladding_k": cladding extinction coefficient,\n
                "cladding_t": cladding thickness (0 = infinite),\n
                "first_n": first mirror refractive index,\n
                "first_k": first mirror extinction coefficient,\n
                "first_t": first mirror thickness (nm),\n
                "second_n": second mirror refractive index,\n
                "second_k": second mirror extinction coefficient,\n
                "second_t": second mirror thickness (nm),\n
                "cavity_n": cavity refractive index,\n
                "cavity_k": cavity extinction coefficient,\n
                "cavity_t": cavity thickness (nm),\n
                "substrate_n": substrate refractive index,\n
                "substrate_k": substrate extinction coefficient,\n
                "substrate_t": substrate thickness (0 = infinite),\n
                "number_of_pairs": number of pairs,\n
                "cavity_index": cavity position,\n
                "symmetrical": "True/False",\n
                "extra_layer": "True/False",\n
                "wavelength_i": wavelength start,\n
                "wavelength_f": wavelength end,\n
                "wavelength_s": wavelength step\n
            }
        Plot settings dictionary containing:
            {
                "width": plot width,\n
                "height": plot height,\n
                "dpi": dots per square inch,\n
                "grid": True/False,\n
                "legend_loc": legend location,\n
                "legend_col": legend column number,\n
                "legend_size": size of legend text,\n
                "axis_fontsize": font size for axis labels,\n
                "label_size": size for tick labels
            }
    file_stem, out_path: string
        File stem for saving. Path to save directory.
    cavity: boolean
        If True, inserts a cavity into the Bragg stack.

    Returns
    -------
    None.

    See Also
    --------
    create_Bragg_stacks
    stack_plots
    Bragg_transfer_matrices
    Bragg_plots

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
    graph_file = Path(f'{out_path}/{file_stem}.png')

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
        extra_layer=simulation_parameters["extra_layer"],
        cavity=cavity)

    plot.stack_plots(
        stack_n=n,
        stack_t=t,
        plot_dict=plot_dict,
        out_path=Path(f'{out_path}/{file_stem}_IndexThickness.png'))

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
    results_path = Path(f'K://Josh/Post_Doc/Python/DBR')
    plot_dict = io.load_json(
        file_path=Path(f'{root}/Plot_Dict.json'))
    sim_params = io.load_json(
        file_path=Path(f'{root}/Simulation_Parameters.json'))

    """ Outputs """
    out_path = Path(
        f'{results_path}/'
        f'n0_{sim_params["cladding_n"]}_'
        f'nc_{sim_params["cavity_n"]}_'
        f'n1_{sim_params["first_n"]}_'
        f'n2_{sim_params["second_n"]}_'
        f'ns_{sim_params["substrate_n"]}_'
        f'MirrorPairs')
    fp.check_dir_exists(dir_path=out_path)

    """ Mirror Pairs """
    mirror_pair_range = range(1, 21, 1)

    for pairs in mirror_pair_range:
        """ Update Simulation Params"""
        sim_params["number_of_pairs"] = pairs

        """ Run Simulation """
        file_stem = f'Bragg_Mirror_Python_MirrorPairs_{pairs}'
        distributed_Bragg_reflector(
            simulation_parameters=sim_params,
            out_path=out_path,
            file_stem=file_stem,
            plot_dict=plot_dict,
            cavity=False)
