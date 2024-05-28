%% Multilayer Bragg Reflector Program

% Benchmarking an older model in an attempt to rewrite this software in
% Python and RCWA alternatives. The code was written in a "add-hoc" manner
% such that there is little-to-no-readability in the original. This program
% serves to fix this problem and then benchmark against the alternatives.

% Clear all previous data
clear all;
%% Initial Parameters

% Here we set the initial parameters for the simulation, this includes the
% refractive indices of the surrounding layers. This includes the cladding
% index, the substrate index, and the effective index of the Bragg Mirror
% pairs. Given that we are going to work layer-by-layer, it is fair to
% assume that the effective index is 0. We can also set the cladding and
% substrate layers to air, if we wish, with an index of 1. We can also se
% the number of mirror pairs and even set the cavity wavelength that we
% want to achieve so that we can get the layer thicknesses. Here, the first
% refractive index (n_first) is that of the first layer of the pair, and
% the second refractive index (n_second) is that of the second layer of the
% pair. Note that the cavity refractive index (n_cavity) can be different.
% We will also set the cavity layer index, this is the layer at which the
% cavity occurs.

n_cladding = 1.0;
n_substrate = 1.0;
n_effective = 0.0;
number_mirror_pairs = 14;
cavity_wavelength = 750;
n_first = 3.6;
n_second = 3.34;
n_cavity = 3.65;
cavity_layer = 7;
wavelength_i = 600;
wavelength_f = 1000;
wavelength_s = 1;
%% Calculations

% Once we set the initial parameters, we can use the standard Bragg Mirror
% equations to calculate the layer thicknesses that we expect to get the
% Bragg response that we want. Note that I have included a commented-out
% variable that allows us to set the thickness that we want too.

t_first = 54;
t_second = 58;
t_cavity = 1000;
%t_first = cavity_wavelength / (4 * n_first);
%t_second = cavity_wavelength / (4 * n_second);
%t_cavity = (cavity_wavelength / 2) / (n_cavity);
%% Set Indices

% Let's set up the index and thickness arrays for the Bragg stack. We will
% loop through the mirror pairs and set the values for x, and x + 1, where
% x is the first of the pair, and x + 1 is the second of the pair.

% Preallocate arrays for refractive indices and thicknesses
n = zeros(1, 2 * number_mirror_pairs);
t = zeros(1, 2 * number_mirror_pairs);

% Set the values for each layer
for x = 1:2:(2 * number_mirror_pairs)
    n(x) = n_first;
    t(x) = t_first;
    if x + 1 <= 2 * number_mirror_pairs
        n(x + 1) = n_second;
        t(x + 1) = t_second;
    end
end

% Set the cavity layer
n(cavity_layer) = n_cavity;
t(cavity_layer) = t_cavity;
%% Loop Wavelengths and Calculate Matrix

% The crux of this simulation is the transfer matrix calculation at each
% layer interface for the entire wavelength range. Log the variables we
% care about so we can see what's going on with the code.

% Number of iterations for wavelength range
number_iterations = (wavelength_f - wavelength_i) / wavelength_s + 1;
matrix_size = [2, 2];

% Preallocate logs
log_matrix = zeros(matrix_size(1), matrix_size(2), number_iterations);
log_k0 = zeros(1, number_iterations);
log_kc = zeros(1, number_iterations);
log_ks = zeros(1, number_iterations);
log_T = zeros(1, number_iterations);
log_R = zeros(1, number_iterations);
log_P = zeros(1, number_iterations);
wavelength = zeros(1, number_iterations);

% Loop over wavelengths
for m = 1:number_iterations
    lambda = wavelength_i + (wavelength_s * (m - 1));
    wavelength(m) = lambda;
    k0 = 2 * pi / lambda;
    log_k0(m) = k0;
    M0 = eye(2);  % Initialize M0 as the identity matrix
    for x = 1:length(n)
        refractive_index = n(x)^2 - n_effective^2;
        kappa = k0 * sqrt(refractive_index);
        kappa_t = kappa * t(x);
        M = [cos(kappa_t), 1i * sin(kappa_t) / kappa; 1i * kappa * sin(kappa_t), cos(kappa_t)];
        M0 = M0 * M;
    end
    log_matrix(:, :, m) = M0;
    kc = k0 * sqrt(n_cladding^2 - n_effective^2);
    log_kc(m) = kc;
    ks = k0 * sqrt(n_substrate^2 - n_effective^2);
    log_ks(m) = ks;
    transmission = 2 * ks / (ks * M0(1, 1) + kc * M0(2, 2) + ks * kc * M0(1, 2) + M0(2, 1));
    reflection = (ks * M0(1, 1) - kc * M0(2, 2) + ks * kc * M0(1, 2) - M0(2, 1)) * transmission / (2 * ks);
    phase = atan(real(transmission) / imag(transmission));
    T = abs(transmission)^2;
    log_T(m) = T;
    R = abs(reflection)^2;
    log_R(m) = R;
    log_P(m) = phase;
end
%% Plot Output

figure('Position', [100, 100, 1200, 1800]);  % [left, bottom, width, height]

% Plot Reflection
subplot(3, 1, 1);
plot(wavelength, log_R, 'LineWidth', 1.5);
title('Reflection Spectrum');
xlim([wavelength_i, wavelength_f]);
ylim([0, 1]);
legend_string = sprintf('t_{first} = %.2f nm, t_{second} = %.2f nm, t_{cavity} = %.2f nm', t_first, t_second, t_cavity);
legend(legend_string, 'Location', 'Best');
xlabel('Wavelength [nm]');
ylabel('Reflection [au]');
grid on;
set(gca, 'FontSize', 12);

% Plot Transmission
subplot(3, 1, 2);
plot(wavelength, log_T, 'LineWidth', 1.5);
title('Transmission Spectrum');
xlim([wavelength_i, wavelength_f]);
ylim([0, 1]);
xlabel('Wavelength [nm]');
ylabel('Transmission [au]');
grid on;
set(gca, 'FontSize', 12);

% Plot Phase
subplot(3, 1, 3);
plot(wavelength, log_P, 'LineWidth', 1.5);
title('Phase Spectrum');
xlim([wavelength_i, wavelength_f]);
xlabel('Wavelength [nm]');
ylabel('Phase [radians]');
grid on;
set(gca, 'FontSize', 12);

% Save the figure
save_path = 'C:\Users\jm1261\Documents\Github\Optical_Simulations\Distributed_Bragg_Reflector\Results\';
filename = sprintf('%sDBR_Cavity_t%.2fnm_t1%.2fnm_t2%.2fnm_MATLAB.png', save_path, t_cavity, t_first, t_second);
%saveas(gcf, filename);
%% Output Logs to Text File
% Specify the text file paths
matrix_log_filename = sprintf('%sDBR_Cavity_t%.2fnm_t1%.2fnm_t2%.2fnm_MATLAB_matrix_log.txt', save_path, t_cavity, t_first, t_second);
other_log_filename = sprintf('%sDBR_Cavity_t%.2fnm_t1%.2fnm_t2%.2fnm_MATLAB_other_log.txt', save_path, t_cavity, t_first, t_second);

% Open the matrix log text file for writing
matrix_log_fileID = fopen(matrix_log_filename, 'w');

% Write the logs of transfer matrices to the matrix log file
for m = 1:number_iterations
    fprintf(matrix_log_fileID, 'Wavelength: %.2f nm\n', wavelength(m));
    fprintf(matrix_log_fileID, 'Matrix:\n');
    fprintf(matrix_log_fileID, '%.6f + %.6fi %.6f + %.6fi\n', real(log_matrix(1,1,m)), imag(log_matrix(1,1,m)), real(log_matrix(1,2,m)), imag(log_matrix(1,2,m)));
    fprintf(matrix_log_fileID, '%.6f + %.6fi %.6f + %.6fi\n\n', real(log_matrix(2,1,m)), imag(log_matrix(2,1,m)), real(log_matrix(2,2,m)), imag(log_matrix(2,2,m)));
end

% Close the matrix log text file
fclose(matrix_log_fileID);

% Open the other log text file for writing
other_log_fileID = fopen(other_log_filename, 'w');

% Write the logs of other parameters to the other log file
fprintf(other_log_fileID, 'Wavelength\tLog of k0\tLog of kc\tLog of ks\tLog of Transmission Coefficients\tLog of Reflection Coefficients\tLog of Phases\n');
for m = 1:number_iterations
    fprintf(other_log_fileID, '%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', wavelength(m), log_k0(m), log_kc(m), log_ks(m), log_T(m), log_R(m), log_P(m));
end

% Close the other log text file
fclose(other_log_fileID);
%% Output Logs

% Display logs for verification
disp('Log of Transfer Matrices:');
disp(log_matrix);
disp('Log of k0:');
disp(log_k0);
disp('Log of kc:');
disp(log_kc);
disp('Log of ks:');
disp(log_ks);
disp('Log of Transmission Coefficients:');
disp(log_T);
disp('Log of Reflection Coefficients:');
disp(log_R);
disp('Log of Phases:');
disp(log_P);