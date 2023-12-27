%% READ TIME SAMPLES AT PROBES PLACED ALONG A WALL_NORMAL LINE

hinfo  =  hdf5info('time_samples.hdf5');

% sampling time
t_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(1));
% wall-normal location of the samples
x_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(5))+1.0;

% sampled velocity components
% each row represents a time instant as dictated by t_smpl
% each column represents a spatial location as dictated by y_smpl
w_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(2));
u_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(3));
v_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(4));


%% Instantaneous velocity plots
figure(1)
hold on
plot(x_smpl, u_smpl(1,:), ':r', 'LineWidth', 2);
plot(x_smpl, v_smpl(1,:), '--b', 'LineWidth', 2);
plot(x_smpl, w_smpl(1,:), '-k', 'LineWidth', 2);
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0, 1]);
legend('$u/u_b$', '$v/u_b$', '$w/u_b$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

%mean value
figure(5)
hold on
plot(x_smpl, mean_x, ':r', 'LineWidth', 2);
plot(x_smpl, mean_y, '--b', 'LineWidth', 2);
plot(x_smpl, mean_z, '-k', 'LineWidth', 2);
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
legend('$u/u_b$', '$v/u_b$', '$w/u_b$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

% Compute the mean and standard deviation of each velocity component at every sampling point
mean_x = mean(u_smpl);
mean_y = mean(v_smpl);
mean_z = mean(w_smpl);
std_x = std(u_smpl);
std_y = std(v_smpl);
std_z = std(w_smpl);

% Plot the results
figure;
subplot(2, 1, 1);
plot(mean_x, 'r');
hold on;
plot(mean_y, 'g');
plot(mean_z, 'b');
title('Velocity Components');
xlabel('Sampling Point');
ylabel('Velocity');
legend('X', 'Y', 'Z');

subplot(2, 1, 2);
plot(std_x, 'r');
hold on;
plot(std_y, 'g');
plot(std_z, 'b');
title('Standard Deviation of Velocity Components');
xlabel('Sampling Point');
ylabel('Standard Deviation');
legend('X', 'Y', 'Z');


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
y0_val = 0.01;    % Displacement height, adjust as needed

% Substitute numerical values into the expressions
U_y_numeric = subs(U_y, [kappa, y0], [kappa_val, y0_val]);
dU_dy_numeric = subs(dU_dy, [kappa, y0], [kappa_val, y0_val]);

% Create function handles for the log-law equation and its gradient
U_y_func = matlabFunction(U_y_numeric);
dU_dy_func = matlabFunction(dU_dy_numeric);

% Create a uniform grid for visualization
y_values = linspace(y0_val, 2, 100); % Adjust the range as needed

% Evaluate the log-law equation and its gradient for numeric values
U_values = U_y_func(y_values);
dU_dy_values = dU_dy_func(y_values);

% Plot both the log-law equation and its gradient
figure;
subplot(2, 1, 1);
plot(y_values, U_values, '-o');
title('Log-Law Equation for Mean Streamwise Velocity');
xlabel('Wall-normal distance (y)');
ylabel('Mean Streamwise Velocity (U)');

subplot(2, 1, 2);
plot(y_values, dU_dy_values, '-o');
title('Analytical Solution for Wall-normal Gradient');
xlabel('Wall-normal distance (y)');
ylabel('Wall-normal Gradient (dU/dy)');

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
