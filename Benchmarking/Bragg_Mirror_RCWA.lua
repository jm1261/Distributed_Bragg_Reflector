-- Create a new simulation, labelled S.
S = S4.NewSimulation()

-- Call the arguments passed to the simulation via the control Python script.
pcall(loadstring(S4.arg))

-- Set the number of harmonics for the Fourier transform.
S:SetNumG(harmonics)

-- Set up cladding and substrate materials.
S:AddMaterial('Cladding', {cladding_n*cladding_n, cladding_k*cladding_k})
S:AddMaterial('Substrate', {substrate_n*substrate_n, substrate_k*substrate_k})

-- Set up Bragg Mirror materials.
S:AddMaterial('First', {first_n*first_n, first_k*first_k})
S:AddMaterial('Second', {second_n*second_n, second_k*second_k})
S:AddMaterial('Cavity', {cavity_n*cavity_n, cavity_k*cavity_k})

-- Split pre- and post- cavity mirror pairs.
pairs_before_cavity = cavity_index - 1

-- Build layers, starting with cladding layer of infinite thickness.
S:AddLayer('TopLayer', 0, 'Cladding')

-- Loop pre-cavity mirror pairs.
for i = 1, pairs_before_cavity do
    S:AddLayer('FirstLayer' .. i, first_t, 'First')
    S:AddLayer('SecondLayer' .. i, second_t, 'Second')
end

-- Add cavity.
S:AddLayer('CavityLayer', cavity_t, 'Cavity')

-- Loop post-cavity mirror pairs.
for i = cavity_index + 1, number_of_pairs - 1 do
    S:AddLayer('SecondLayer' .. i, second_t, 'Second')
    S:AddLayer('FirstLayer' .. i, first_t, 'First')
end

-- Add the final mirror layer
S:AddLayer('SecondLayer', second_t, 'Second')

-- Add the substrate layer at the bottom
S:AddLayer('BottomLayer', 0, 'Substrate')

-- Set the light source parameters. {phi[0, 180], theta[0, 360]}{s-amp, s-phase}
--{p-amp, p-phase}.
S:SetExcitationPlanewave(
    {0, 0},
    {1, 0},
    {0, 0})

-- Loop wavelength range and get transmission, reflection, and phase.
for lambda = wavelength_i, wavelength_f, wavelength_s do
    freq = 1/lambda
    S:SetFrequency(freq)

    -- Get transmitted power flux
    transmission = S:GetPowerFlux('BottomLayer')
    
    -- Get incident and reflected power flux
    inc, reflection = S:GetPowerFlux('TopLayer', 10)

    print(lambda, transmission, -reflection)
end