%% Demo: explicit model, scalar model parameter, scalar measurement point
% Different use cases of the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> are illustrated here. This application
% example is an explicit model with one model parameter and one-dimensional
% measurement points.


%% Create the model object
x = sym('x');                               % The independent variable
p = sym('p');                               % The model parameter
f = p * x;                                  % The model function
model = model_explicit(f, p, x);            % Create the model object using model_explicit

%% Create the solver object
p = 1                                       % True parameter of the model
p0 = (1 + rand()) * p                       % Guessed parameter value

n = 5;                                      % Number of different selectable measurements
x_var = (0:1/(n-1):1)'                      % Selectable measurements
v_var = 10^-2 * ones(1, n)'                 % Variances of measurement results at these measurements

sol = solver(model, p0, x_var, v_var);      % Create the solver object


%% Calculate optimal measurements
max = 3;                                    % Maximal number of measurements to choose
x_opt = sol.get_optimal_measurements(max)   % Calculate the optimal measurements of the selectable measurements


%% Calculate quality of measurements
% The smaller the value, the better the quality.
w_opt = sol.get_optimal_weights(max)            % Calculate the optimal weights of the selectable measurements
quality_opt = sol.get_quality(w_opt)            % Calculate quality resulting from optimal measurements
w_subopt = [ones(max, 1); zeros(n-max, 1)]      % Suboptimal weights
quality_subopt = sol.get_quality(w_subopt)      % Calculate quality resulting from suboptimal measurements


%% Estimate model parameter from accomplished measurements
m = 5;                                                                          % Number of accomplished measurements
x_fix = x_opt;                                                                  % Accomplished measurements
v_fix = v_var(w_opt);                                                           % Variances of measurement results at these measurements
eta = model_util.get_fictitious_measurement_results(model, p, x_fix, v_fix);    % Measurement results of the accomplished measurements
sol.set_accomplished_measurements(x_fix, v_fix, eta);                           % Pass accomplished measurements to the solver object
p_lb = 0;                                                                       % Lower bound of model parameter
p_ub = 2;                                                                       % Upper bound of model parameter
p_opt = sol.get_optimal_parameters(p_lb, p_ub)                                  % Optimize model parameter from accomplished measurements


%% Calculate gain of additional measurements
sol.set_initial_parameter_estimation(p_opt);   % Update parameter estimation
w_opt = sol.get_optimal_weights(max)';         % Calculate the optimal weights of the selectable measurements
quality_old = sol.get_quality(zeros(n, 1))     % Calculate quality without additional measurements
quality_new = sol.get_quality(w_opt)           % Calculate quality resulting from optimal additional measurements


%% Calculate optimal measurements with constraints
% We are constraining the choice of measurements in such a way that
% distance between two chosen measurements has to be at least 0.75.
A = diag(ones(n, 1)) + diag(ones(n-1, 1), 1) + diag(ones(n-2, 1), 2)    % Matrix for the constraints of the measurements
b = ones(n, 1)                                                          % Vector for the constraints of the measurements
x_opt = sol.get_optimal_measurements(A, b)                              % Calculate the optimal measurements of the selectable measurements considering the constraints

