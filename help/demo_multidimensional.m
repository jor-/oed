%% Multidimensional demo
% Different use cases of the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> are illustrated here with an example
% with an multidimensional model parameter vector and multidimensional measurements.


%% Create the model and the solver object
p = [0, 1];                                         % True parameters of the model
p0 = p + rand(size(p)) - 0.5                        % Guessed parameter values

n_x = 3;                                            % Number of different selectable measurements for the x variable
n_y = 4;                                            % Number of different selectable measurements for the y variable
x_var = (0:1/(n_x-1):1);                            % Selectable measurements for the x variable
y_var = (0:1/(n_y-1):1);                            % Selectable measurements for the y variable
[x_var_tmp, y_var_tmp] = meshgrid(x_var, y_var);    % Temporarily variables for combination of both selectable measurements
t_var = [x_var_tmp(:) y_var_tmp(:)]                 % Selectable measurements for both variables
n = n_x * n_y;                                      % Number of different selectable measurements for both variables
v_var = 10^-2 * ones(length(t_var), 1)              % Variances of measurement results at these measurements

model = model_explicit('a*x^2 + b*y', {'a', 'b'}, {'x', 'y'});  % Create the model object
sol = solver(model, p0, t_var, v_var);                          % Create the solver object


%% Calculate optimal measurements
max = 3;                                    % Maximal number of measurements to choose
t_opt = sol.get_optimal_measurements(max)   % Calculate the optimal measurements of the selectable measurements


%% Calculate quality of measurements
% The smaller the value, the better the quality.
w_opt = sol.get_optimal_weights(max)            % Calculate the optimal weights of the selectable measurements
quality_opt = sol.get_quality(w_opt)            % Calculate quality resulting from optimal measurements
w_subopt = [ones(max, 1); zeros(n-max, 1)]      % Suboptimal weights
quality_subopt = sol.get_quality(w_subopt)      % Calculate quality resulting from suboptimal measurements


%% Estimate model parameters from accomplished measurements
m = 5;                                                                          % Number of accomplished measurements
t_fix = t_opt;                                                                  % Accomplished measurements
v_fix = v_var(w_opt);                                                           % Variances of measurement results at these measurements
eta = model_util.get_fictitious_measurement_results(model, p, t_fix, v_fix);    % Measurement results of the accomplished measurements
sol.set_accomplished_measurements(t_fix, v_fix, eta);                           % Pass accomplished measurements to the solver object
p_lb = [-1, 0];                                                                 % Lower bounds of model parameters
p_ub = [1, 2];                                                                  % Upper bounds of model parameters
p_opt = sol.get_optimal_parameters(p_lb, p_ub)                                  % Optimize model parameter from accomplished measurements


%% Calculate gain of additional measurements
sol.set_initial_parameter_estimation(p_opt);   % Update parameter estimation
w_opt = sol.get_optimal_weights(max);          % Calculate the optimal weights of the selectable measurements
quality_old = sol.get_quality(zeros(n, 1))     % Calculate quality without additional measurements
quality_new = sol.get_quality(w_opt)           % Calculate quality resulting from optimal additional measurements


%% Calculate optimal measurements with constraints
% We are constraining the choice of measurements in such a way that
% different measurements should have different values in the x variable.
A_tmp = diag(ones(4, 1)) + diag(ones(3, 1), 1) + diag(ones(2, 1), 2) + diag(ones(1, 1), 3)      % Temporarily matrix for the constraints of the measurements
A = blkdiag(A_tmp, A_tmp, A_tmp)            % Matrix for the constraints of the measurements
b = ones(n, 1)                              % Vector for the constraints of the measurements
t_opt = sol.get_optimal_measurements(A, b)  % Calculate the optimal measurements of the selectable measurements considering the constraints
