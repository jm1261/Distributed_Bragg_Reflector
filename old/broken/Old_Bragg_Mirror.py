import numpy as np
import matplotlib.pyplot as plt
import cmath

## Define outer layer of refractive index ##
n_clad = 1.0 ## cladding ##
n_sub = 1.0 ## substrate ##
n_eff = 0.0 ## effective ##

## Define parameters of repeating structure ##
n1 = 3.11
t1 = 52
n2 = 3.53
t2 = 46

## Set cavity gain ##
gain = 0

## Define layers of the material ##
n = []
t = []
for a in range(0, 80, 2):
    n.append(n1)
    t.append(t1)
    n.append(n2)
    t.append(t2)

## Now set the defect ##
n[15] = 3.53 + 0j
t[15] = 92

wavs = []
reflection = []
transmission = []
phase_array = []
for b in range(2000):
    wavelength = 550 + 0.1*b
    wavs.append(wavelength)
    k0 = (2 * cmath.pi) / wavelength
    M0 = [[1, 0],
          [0,1]]

    for c in range(len(n)):
        ref_i = (n[c]*n[c]) - (n_eff*n_eff)
        nk = k0 * cmath.sqrt(ref_i)
        nkt = nk * t[c]
        M = [[cmath.cos(nkt), 1j/nk * cmath.sin(nkt)],
             [1j*nk * cmath.sin(nkt), cmath.cos(nkt)]]
        M1 = np.matmul(M0, M)
        M0 = M1

    k_clad = k0 * cmath.sqrt(n_clad * n_clad
                             - n_eff * n_eff)
    k_sub = k0 * cmath.sqrt(n_sub * n_sub
                            - n_eff * n_eff)

    ## Transmission ##
    trans = (
    ((2 * k_sub) / (k_sub * M0[0][0])) +
    (k_clad * M0[1][1]) +
    (k_sub * k_clad * M0[0][1]) +
    (M0[1][0])
    )
    ## Reflection ##
    rflct = (
    (k_sub * M0[0][0]) -
    (k_clad * M0[1][1]) +
    (k_sub * k_clad * M0[0][1]) -
    (M0[1][0] * (trans / (2 * k_clad)))
    )

    ## Phase ##
    phase = cmath.atan(trans.real / trans.imag)

    ## Update the matrices##
    transmission.append(trans * trans.conjugate())
    reflection.append(rflct * rflct.conjugate())
    phase_array.append(phase)


plt.plot(wavs, reflection)
plt.show()
plt.clf()

plt.plot(wavs, transmission)
plt.show()
plt.clf()
