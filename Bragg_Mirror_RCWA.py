# %%
import numpy as np
import src.dataIO as io
import src.analysis as anal
import src.plotting as plot

from pathlib import Path

# %%
simulation_parameters = {
    "harmonics": 15,
    "cladding_n": 1,
    "cladding_k": 0,
    "cladding_t": 0,
    "first_n": 3.6,
    "first_k": 0,
    "first_t": 54,
    "second_n": 3.34,
    "second_k": 0,
    "second_t": 58,
    "cavity_n": 3.65,
    "cavity_k": 0,
    "cavity_t": 1000,
    "substrate_n": 1,
    "substrate_k": 0,
    "substrate_t": 0,
    "number_of_pairs": 14,
    "cavity_index": 5,
    "wavelength_i": 600,
    "wavelength_f": 1000,
    "wavelength_s": 1}

# %%
root = Path(f'C:/Users/jm1261/Documents/Github/Optical_Simulations/Distributed_Bragg_Reflector')
results_path = Path(f'{root}/Results')
plot_dict = io.load_json(file_path=Path(f'{root}/Plot_Dict.json'))

# %%
filename = (
    f'DBR_Cavity_t{simulation_parameters["cavity_t"]}nm'
    f'_t1{simulation_parameters["first_t"]}nm'
    f'_t2{simulation_parameters["second_t"]}nm')
log_file = Path(f'{results_path}/{filename}_RCWA_log.txt')
graph_file = Path(f'{results_path}/{filename}_RCWA.png')

# %%
arg_string = anal.argument_string(
    names=[key for key, _ in simulation_parameters.items()],
    values=[value for _, value in simulation_parameters.items()])
print(arg_string)

# %%
S4_output = anal.S4_RCWA(
    lua_script=Path(f'{root}/Bragg_Mirror_RCWA.lua'),
    argument_string=arg_string)
wavelength, transmission, reflection = anal.read_S4_output(
    process_string=S4_output)

# %%
data = [
    [wav, transmission[i], reflection[i]]
    for i, wav in enumerate(wavelength)]
np.savetxt(
    log_file,
    data,
    delimiter=',')

# %%
plot.S4_plot(
    wavelength=wavelength,
    transmission=transmission,
    reflection=reflection,
    out_path=graph_file,
    plot_dict=plot_dict)

# %%
plot.Bragg_plots(
    wavelength=wavelength,
    transmission=transmission,
    reflection=reflection,
    phase=np.zeros(len(wavelength)),
    out_path=graph_file,
    plot_dict=plot_dict)
