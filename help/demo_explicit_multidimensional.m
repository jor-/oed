%% Demo: explicit model, multidimensional model parameter vector, multidimensional measurement points
% Different use cases of the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> are illustrated here. This application
% example is an explicit model with multidimensional model parameter vector
% and multidimensional measurement points.


%% Create the model object
t = sym('t');                                       % Create symbolic variable for model
s = sym('s');                                       % Create symbolic variable for model
a = sym('a');                                       % Create symbolic variable for model
b = sym('b');                                       % Create symbolic variable for model
x = [t, s];                                         % The independent variables
p = [a, b];                                         % The model parameters
f = a*t^2 + b*s;                                    % The model function
model = model_explicit(f, p, x);                    % Create the model object using model_explicit

%% Create the solver object
p = [0; 1]                                          % True parameters of the model
p0 = p + rand(size(p)) - 0.5                        % Guessed parameter values

n_t = 3;                                            % Number of different selectable measurements for the t variable
n_s = 4;                                            % Number of different selectable measurements for the s variable
t_var = (0:1/(n_t-1):1);                            % Selectable measurements for the x variable
s_var = (0:1/(n_s-1):1);                            % Selectable measurements for the y variable
[t_var_tmp, s_var_tmp] = meshgrid(t_var, s_var);    % Temporarily variables for combination of both selectable measurements
x_var = [t_var_tmp(:) s_var_tmp(:)]                 % Selectable measurements for both variables
n = n_t * n_s;                                      % Number of different selectable measurements for both variables
v_var = 10^-2 * ones(length(x_var), 1)              % Variances of measurement results at these measurements

sol = solver(model, p0, x_var, v_var);              % Create the solver object


%% Calculate optimal measurements
max = 3;                                    % Maximal number of measurements to choose
x_opt = sol.get_optimal_measurements(max)   % Calculate the optimal measurements of the selectable measurements


%% Calculate quality of measurements
% The smaller the value, the better the quality.
w_opt = sol.get_optimal_weights(max)            % Calculate the optimal weights of the selectable measurements
quality_opt = sol.get_quality(w_opt)            % Calculate quality resulting from optimal measurements
w_subopt = [ones(max, 1); zeros(n-max, 1)]      % Suboptimal weights
quality_subopt = sol.get_quality(w_subopt)      % Calculate quality resulting from suboptimal measurements


%% Estimate model parameters from accomplished measurements
m = 5;                                                                          % Number of accomplished measurements
x_fix = x_opt;                                                                  % Accomplished measurements
v_fix = v_var(w_opt);                                                           % Variances of measurement results at these measurements
eta = model_util.get_fictitious_measurement_results(model, p, x_fix, v_fix);    % Measurement results of the accomplished measurements
sol.set_accomplished_measurements(x_fix, v_fix, eta);                           % Pass accomplished measurements to the solver object
p_lb = [-1; 0];                                                                 % Lower bounds of model parameters
p_ub = [1; 2];                                                                  % Upper bounds of model parameters
p_opt = sol.get_optimal_parameters(p_lb, p_ub)                                  % Optimize model parameter from accomplished measurements


%% Calculate gain of additional measurements
sol.set_initial_parameter_estimation(p_opt);   % Update parameter estimation
w_opt = sol.get_optimal_weights(max);          % Calculate the optimal weights of the selectable measurements
quality_old = sol.get_quality(zeros(n, 1))     % Calculate quality without additional measurements
quality_new = sol.get_quality(w_opt)           % Calculate quality resulting from optimal additional measurements


%% Calculate optimal measurements with constraints
% We are constraining the choice of measurements in such a way that different
% measurements should have different values in the first independet variable (t).
A_tmp = diag(ones(4, 1)) + diag(ones(3, 1), 1) + diag(ones(2, 1), 2) + diag(ones(1, 1), 3)      % Temporarily matrix for the constraints of the measurements
A = blkdiag(A_tmp, A_tmp, A_tmp)            % Matrix for the constraints of the measurements
b = ones(n, 1)                              % Vector for the constraints of the measurements
x_opt = sol.get_optimal_measurements(A, b)  % Calculate the optimal measurements of the selectable measurements considering the constraints

