# %%
import numpy as np
import src.dataIO as io
import src.plotting as plot

from pathlib import Path

#%%
sim_params = {
    "n_effective": 0,
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
    "cavity_index": 6,
    "wavelength_i": 600,
    "wavelength_f": 1000,
    "wavelength_s": 0.1}

# %%
root = Path(f'C:/Users/jm1261/Documents/Github/Optical_Simulations/Distributed_Bragg_Reflector')
results_path = Path(f'{root}/Results')
plot_dict = io.load_json(file_path=Path(f'{root}/Plot_Dict.json'))

# %%
filename = (
    f'DBR_Cavity_t{sim_params["cavity_t"]}nm'
    f'_t1{sim_params["first_t"]}nm'
    f'_t2{sim_params["second_t"]}nm')
log_file = Path(f'{results_path}/{filename}_Python_log.txt')
graph_file = Path(f'{results_path}/{filename}_Python.png')

# %%
transmission = []
reflection = []
phase = []

# %%
n = (
    [sim_params["first_n"], sim_params["second_n"]] *
    sim_params["number_of_pairs"])
t = (
    [sim_params["first_t"], sim_params["second_t"]] *
    sim_params["number_of_pairs"])

n[sim_params["cavity_index"]] = sim_params["cavity_n"]
t[sim_params["cavity_index"]] = sim_params["cavity_t"]

# %%
wavelength_range = np.arange(
    sim_params["wavelength_i"],
    sim_params["wavelength_f"] + sim_params["wavelength_s"],
    sim_params["wavelength_s"])

# %%
for wav in wavelength_range:
    k0 = (2 * np.pi) / wav
    M0 = [[1, 0], [0, 1]]

    for i in range(sim_params["number_of_pairs"] * 2):
        refractive_index = (
            (n[i] ** 2) - (sim_params["n_effective"] ** 2))
        nk = np.sqrt(refractive_index) * k0
        nkt = nk * t[i]

        M = [
            [np.cos(nkt), (1j / nk) * np.sin(nkt)],
            [1j * nk * np.sin(nkt), np.cos(nkt)]]
        M0 = np.matmul(M0, M)

    kc = k0 * np.sqrt(
        (sim_params["cladding_n"] ** 2) -
        (sim_params["n_effective"] ** 2))
    ks = k0 * np.sqrt(
        (sim_params["substrate_n"] ** 2) -
        (sim_params["n_effective"] ** 2))

    T = (
        (2 * ks) /
        (
            (ks * M0[0][0]) +
            (kc * M0[1][1]) +
            (ks * kc * M0[0][1]) +
            (M0[1][0])))
    transmission.append(T * T.conjugate())

    R = (
        (
            (ks * M0[0][0]) -
            (kc * M0[1][1]) +
            (ks * kc * M0[0][1]) -
            (M0[1][0])) *
        (T / (2 * ks)))
    reflection.append(R * R.conjugate())

    P = np.arctan(T.real / T.imag)
    phase.append(P)

# %%
data = [
    [wav, transmission[i], reflection[i], phase[i]]
    for i, wav in enumerate(wavelength_range)]
np.savetxt(
    log_file,
    data,
    delimiter=',')

# %%
plot.Bragg_plots(
    wavelength=wavelength_range,
    transmission=transmission,
    reflection=reflection,
    phase=phase,
    out_path=graph_file,
    plot_dict=plot_dict)

# %%
