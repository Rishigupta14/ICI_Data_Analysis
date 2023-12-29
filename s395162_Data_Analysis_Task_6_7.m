clear all
close all

FileName = 'time_samples.hdf5';

File = h5info(FileName);

% Sampling time
t_smpl = h5read(File.Filename, '/t');

% Wall-normal location of the samples
x_smpl = h5read(File.Filename, '/y');

% Sampled velocity components
w_smpl = h5read(File.Filename, '/w'); % Streamwise
u_smpl = h5read(File.Filename, '/u'); % Wall Normal
v_smpl = h5read(File.Filename, '/v'); % Spanwise

% Find the index corresponding to x=0.06
x_index = find(x_smpl <= 0.06, 1, "last");

W_x_p06 = interp1(x_smpl, w_smpl.', 0.06, 'linear', 'extrap').';
W_x_n06 = interp1(x_smpl, w_smpl.', -0.06, 'linear', 'extrap').';
figure                                                          % Explore Data
surfc(x_smpl, t_smpl, w_smpl, 'EdgeColor', 'interp')
hold on
plot3(0.06 * ones(size(t_smpl)), t_smpl, W_x_p06, '-r', 'LineWidth', 1)
plot3(-0.06 * ones(size(t_smpl)), t_smpl, W_x_n06, '-g', 'LineWidth', 3)
hold off
grid
colormap(turbo)
xlabel('x\_smpl')
ylabel('t\_smpl')
zlabel('w\_smpl')
view(60, 30)

Ts = mean(diff(t_smpl));
Fs = 1/Ts;
Fn = Fs/2;
L = numel(t_smpl); NFFT = 2^nextpow2(L);
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

% Extract wall-normal velocity signal at x=0.06
u_signal = u_smpl(:, x_index);

% Ensure the vectors have the same length
min_length = min(length(t_smpl), length(u_signal));
t_smpl = t_smpl(1:min_length);
u_signal = u_signal(1:min_length);

% Perform FFT analysis
delta_t = mean(diff(t_smpl)); % Average time step
Fs = 1 / delta_t; % Sampling frequency
L = length(u_signal);
frequencies = Fs * (0:(L-1)) / L; % Corrected the frequency calculation
fft_result = fft(u_signal);

% Ensure fft_result has a valid length
valid_length = min(L, floor(L/2));
fft_amplitude = 2 / L * abs(fft_result(1:valid_length));

% Plot default FFT spectrum
figure;
subplot(2, 1, 1);
plot(t_smpl, u_signal);
title('Wall-normal Velocity Signal at x=0.06');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

subplot(2, 1, 2);
plot(Fv, abs(FTw(Iv))*2, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum - Default Plot');
grid on;
hold off;

% Logarithmic plot
figure;
semilogx(Fv, abs(FTw(Iv))*2, 'LineWidth', 2);
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

    % Apply window function
    hann_window = 0.5 * (1 - cos(2 * pi * (1:L) / (L + 1)));
    w_windowed = u_signal .* hann_window.';

    % Perform FFT analysis on windowed signal
    fft_result_windowed = fft(w_windowed);
    fft_amplitude_windowed = 2/L * abs(fft_result_windowed(1:valid_length));

    % Truncate Fv to match the size of fft_amplitude_windowed
    Fv_truncated = Fv(1:numel(fft_amplitude_windowed));

    % Plot FFT spectrum with Hann window
    figure;
    % Plot the first numel(fft_amplitude_windowed) elements of Fv and fft_amplitude_windowed
    subplot(2, 1, 1);
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
u_smoothed = movmean(u_signal, window_size);

% Plot original and smoothed signals
figure;
subplot(2, 1, 1);
plot(t_smpl, u_signal, 'LineWidth', 2, 'DisplayName', 'Original Signal');
title('Original and Smoothed Signals');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('show');
grid on;
hold off;

subplot(2, 1, 2);
plot(t_smpl, u_smoothed, 'r--', 'LineWidth', 2, 'DisplayName', 'Smoothed Signal');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('show');
grid on;
hold off;

% Extract the wall-normal velocity signal at the x=0.06 location
signal = u_signal;

% Determine if the signal has a normal distribution and estimate the parameters
figure;
histfit(signal);
pd = fitdist(signal, 'normal');
title(['Histogram of the signal with a normal density fit. Mean: ', num2str(pd.mu), ', Standard deviation: ', num2str(pd.sigma)]);

% Evaluate the statistical significance of the correlation between the two velocity components
corr_coeff = corrcoef(u_signal);
p_value_ks = kstest2(u_signal, u_signal);
p_value_lillie = lillietest(u_signal);
disp(['The correlation coefficient between the two velocity components is: ', num2str(corr_coeff)]);
disp(['The p-value of the Kolmogorov-Smirnov test is: ' num2str(p_value_ks)]);
disp(['The p-value of the Lilliefors test is: ' num2str(p_value_lillie)]);
