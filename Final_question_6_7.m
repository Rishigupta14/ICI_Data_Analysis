
% Symbolic variables
syms U y kappa y0

% Log-law equation
U_y = (1/kappa) * log(y / y0);

% Compute the gradient
dU_dy = diff(U_y, y);

% Display the symbolic expression for the gradient
disp('Analytical expression for the wall-normal gradient:')
disp(dU_dy);

% Specify numerical values for the parameters
kappa_val = 0.41; % von Kármán constant
y0_val = 15;      % Displacement height, adjust as needed

% Substitute numerical values into the expressions
U_y_numeric = subs(U_y, [kappa, y0], [kappa_val, y0_val]);
dU_dy_numeric = subs(dU_dy, [kappa, y0], [kappa_val, y0_val]);

% Create function handles for the log-law equation and its gradient
U_y_func = matlabFunction(U_y_numeric);
dU_dy_func = matlabFunction(dU_dy_numeric);

% Grid convergence analysis
num_grids = 5; % Number of grids for analysis
errors = zeros(num_grids, 1);
grid_sizes = zeros(num_grids, 1);

for i = 1:num_grids
    % Create a uniform grid for visualization
    grid_size = 10 * 2^i; % Adjust the initial grid size and refinement factor as needed
    y_values = linspace(y0_val, 2, grid_size);
    
    % Evaluate the log-law equation and its gradient for numeric values
    U_values = U_y_func(y_values);
    dU_dy_values = dU_dy_func(y_values);
    
    % Compute numerical approximation of the gradient (use finite differences)
    dU_dy_approx = diff(U_values) ./ diff(y_values);
    
    % Compute error
    errors(i) = norm(dU_dy_values(1:end-1) - dU_dy_approx);
    
    % Save grid size for Richardson extrapolation
    grid_sizes(i) = grid_size;
end

% Richardson extrapolation
p = polyfit(log(1./grid_sizes), log(errors), 1);
extrapolated_error = exp(polyval(p, 0));

% Plot results
figure;
loglog(1./grid_sizes, errors, '-o');
hold on;
loglog(1./grid_sizes(end), extrapolated_error, 'rx', 'MarkerSize', 10);
title('Grid Convergence Analysis');
xlabel('1/Grid Size');
ylabel('Error');
legend('Error vs. 1/Grid Size', 'Richardson Extrapolation');
grid on;

disp('Extrapolated error using Richardson extrapolation:')
disp(extrapolated_error);
clear all
close all

FileName = 'time_samples.hdf5';

File = h5info(FileName);

% Sampling time

% t_smpl = h5read(h5info().GroupHierarchy.Datasets.(1))
t_smpl = h5read(File.Filename,'/t');

% Wall-normal location of the samples
% x_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(5)) + 1.0;
x_smpl = h5read(File.Filename,'/y');

% Sampled velocity components
% w_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(2));
w_smpl = h5read(File.Filename,'/w');

% u_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(3));
u_smpl = h5read(File.Filename,'/u');

% v_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(4));
v_smpl = h5read(File.Filename,'/v');

% Find the index corresponding to x=0.06
x_index = find(x_smpl <= 0.06, 1, "last");

W_x_p06 = interp1(x_smpl, w_smpl.', 0.06, 'linear', 'extrap').';
W_x_n06 = interp1(x_smpl, w_smpl.', -0.06, 'linear', 'extrap').';
figure                                                          % Explore Data
surfc(x_smpl, t_smpl, w_smpl, 'EdgeColor','interp')
hold on
plot3(0.06*ones(size(t_smpl)), t_smpl, W_x_p06, '-r', 'LineWidth',1)
plot3(-0.06*ones(size(t_smpl)), t_smpl, W_x_n06, '-g', 'LineWidth',3)
hold off
grid
colormap(turbo)
xlabel('x\_smpl')
ylabel('t\_smpl')
zlabel('w\_smpl')
view(60, 30)

% t_stats = [mean(diff(t_smpl))  std(diff(t_smpl))]

Ts = mean(diff(t_smpl));
Fs = 1/Ts;
Fn = Fs/2;
L = numel(t_smpl);NFFT = 2^nextpow2(L);
FTw = fft((W_x_n06 - mean(W_x_n06)).*hann(L), NFFT)/L;
Fv = linspace(0, 1, NFFT/2+1)*Fn;
Iv = 1:numel(Fv);

[FTw_max, idx] = max(abs(FTw(Iv))*2);

figure
plot(Fv, abs(FTw(Iv))*2)
grid
xlabel('Frequency')
ylabel('Magnitude')
title('Fourier Transform')
xlim([0 6])
text(Fv(idx), FTw_max, sprintf('\\leftarrow Magnitude: %.6f\n     Frequency: %.6f',FTw_max,Fv(idx)), 'Vert','top')

%fprintf('\n\nStopped at line 73.  Remove this line and the following %return’ to continue.\n\n')
%return

% Extract wall-normal velocity signal at x=0.06
w_signal = w_smpl(:, x_index);

% Ensure the vectors have the same length
min_length = min(length(t_smpl), length(w_signal));
t_smpl = t_smpl(1:min_length);
w_signal = w_signal(1:min_length);

% Perform FFT analysis
delta_t = mean(diff(t_smpl)); % Average time step
Fs = 1 / delta_t; % Sampling frequency
L = length(w_signal);
frequencies = Fs * (0:(L-1)) / L; % Corrected the frequency calculation
fft_result = fft(w_signal);

% Ensure fft_result has valid length
valid_length = min(L, floor(L/2));
fft_amplitude = 2 / L * abs(fft_result(1:valid_length));

% Plot default FFT spectrum
figure;
subplot(2, 1, 1);
plot(t_smpl, w_signal);
title('Wall-normal Velocity Signal at x=0.06');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

subplot(2, 1, 2);
plot(Fv, abs(FTw(Iv))*2, 'LineWidth', 2);
%plot(frequencies, fft_amplitude, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum - Default Plot');
grid on;
hold off;

% Logarithmic plot
figure;
semilogx(Fv, abs(FTw(Iv))*2, 'LineWidth', 2);
%semilogx(frequencies, fft_amplitude, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum - Logarithmic Plot');
grid on;
hold off;

% Apply window functions and plot FFT spectra
window_functions = {'rectwin', 'hamming', 'hann', 'blackman'};
figure;
for i = 1:length(window_functions)
    window_function = window_functions{i};
     window = str2func(window_function);

% Apply  window function
hann_window = 0.5 * (1 - cos(2 * pi * (1:L) / (L + 1)));
w_windowed = w_signal .* hann_window.';

% Perform FFT analysis on windowed signal
fft_result_windowed = fft(w_windowed);
fft_amplitude_windowed = 2/L * abs(fft_result_windowed(1:valid_length));

% Truncate Fv to match the size of fft_amplitude_windowed
Fv_truncated = Fv(1:numel(fft_amplitude_windowed));

% Plot FFT spectrum with Hann window
figure;
% Plot the first numel(fft_amplitude_windowed) elements of Fv and fft_amplitude_windowed
plot(Fv_truncated, fft_amplitude_windowed, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum with Hann Window');
grid on;
hold off;
end



xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum with Different Window Functions');
legend('show');
grid on;
hold off;

% Moving average filtering
window_size = 10; % Adjust window size as needed
w_smoothed = movmean(w_signal, window_size);

% Plot original and smoothed signals
figure;
subplot(2, 1, 1);
plot(t_smpl, w_signal, 'LineWidth', 2, 'DisplayName', 'Original Signal');
title('Original and Smoothed Signals');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('show');
grid on;
hold off;

subplot(2, 1, 2);
plot(t_smpl, w_smoothed, 'r--', 'LineWidth', 2, 'DisplayName', 'Smoothed Signal');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('show');
grid on;
hold off;

% Extract the wall-normal velocity signal at the x=0.06 location
signal = w_signal;

% Determine if the signal has a normal distribution and estimate the parameters
figure;
histfit(signal);
pd = fitdist(signal, 'normal');
title(['Histogram of the signal with a normal density fit. Mean: ', num2str(pd.mu), ', Standard deviation: ', num2str(pd.sigma)]);

% Evaluate the statistical significance of the correlation between the two velocity components
corr_coeff = corrcoef(w_signal);
p_value_ks = kstest2(w_signal, w_signal);
p_value_lillie = lillietest(w_signal);
disp(['The correlation coefficient between the two velocity components is: ', num2str(corr_coeff)]);
disp(['The p-value of the Kolmogorov-Smirnov test is: ' num2str(p_value_ks)]);
disp(['The p-value of the Lilliefors test is: ' num2str(p_value_lillie)]);