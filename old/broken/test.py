# %%
import numpy as np
import matplotlib.pyplot as plt
# %%
n_cladding = 1.0
n_substrate = 1.0
n_effective = 0.0
# %%
wavelength_range = range(700, 850, 1)
# %%
fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=[15, 9],
    tight_layout=True)
# %%
T = []
R = []
P = []
resonance_location = []
# %%
number_mirror_pairs = 28
n_high = 3.6
t_high = 54
n_low = 3.34
t_low = 58
cavity_n = 3.65
cavity_t = 100
# %%
n = [n_high, n_low] * number_mirror_pairs
t = [t_high, t_low] * number_mirror_pairs
n[6] = cavity_n
t[6] = cavity_t
# %%
T.clear()
R.clear()
P.clear()
# %%
for wav in wavelength_range:
    k0 = (2 * np.pi) / wav
    M0 = [[1, 0], [0, 1]]

    for i in range(number_mirror_pairs * 2):
        refractive_index = (n[i] ** 2) - (n_effective ** 2)
        kappa = k0 * np.sqrt(refractive_index)
        arg = kappa * t[i]
        M = [
            [np.cos(arg),                   ((1j * np.sin(arg)) / kappa)],
            [(1j * kappa * np.sin(arg)),    np.cos(arg)]]
        M0 = np.matmul(M0, M)

    kc = k0 * np.sqrt((n_cladding ** 2) - (n_effective ** 2))
    ks = k0 * np.sqrt((n_substrate ** 2) - (n_effective ** 2))

    transmission = (
        (2 * ks) /
        (
            (ks * M0[0][0]) +
            (kc * M0[1][1]) +
            (ks * kc * M0[0][1]) +
            (M0[1][0])))

    reflection = (
        (
            (ks * M0[0][0]) -
            (kc * M0[1][1]) +
            (ks * kc * M0[0][1]) -
            (M0[1][0])) *
        (transmission / (2 * ks)))

    phase = np.arctan(transmission.real / transmission.imag)

    T.append(transmission * transmission.conjugate())
    R.append(reflection * reflection.conjugate())
    P.append(phase)
# %%
y = [1 - r.real for r in R]
resonance_location.append(wavelength_range[np.argmax(y)])
# %%
ax.plot(
    wavelength_range,
    y,
    linewidth=4)
ax.set_ylim(0, 1)
# %%
plt.show()
print(M0)