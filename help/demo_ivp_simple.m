%% Demo: inital value problem model (multidimensional model parameter vector, scalar measurement points)
% Different use cases of the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> are illustrated here. The applicaion
% example is an inital value problem as model with multidimensional model
% parameter vector and multidimensional measurements.


%% Create the model object
t = sym('t');                                   % The independent variable
y = sym('y');                                   % The dependent variable
a = sym('a');                                   % A model parameter
y0 = sym('y0');                                 % The inital value as another model parameter (this could also be a fix scalar value)
p = [a; y0];                                    % All model parameters as a vector
f = y + a * t;                                  % The differential equation of the inital value problem
t_interval = [0; 10];                           % The first value is the initial (time) point (y(t_span(1)) = y0), the second value is an upper bound up to where the model is maximally evaluated 
model = model_ivp(f, p, y, y0, t, t_interval);  % Create the model object using model_ivp


%% Create the solver object
p = [2; 1]                                  % True parameters of the model
p0 = p + rand(size(p)) - 0.5                % Guessed parameter values

n = 5;                                      % Number of different selectable measurements
t_var = (0:10/(n-1):10)'                    % Selectable measurements
v_var = 10^-2 * ones(1, n)'                 % Variances of measurement results at these measurements

sol = solver(model, p0, t_var, v_var);      % Create the solver object


%% Calculate optimal measurements
max = 2;                                    % Maximal number of measurements to choose
t_opt = sol.get_optimal_measurements(max)   % Calculate the optimal measurements of the selectable measurements


%% Calculate quality of measurements
% The smaller the value, the better the quality.
sol.set_option('edo_algorithm', 'direct');      % Use direct solver in calculation of optimal weights
w_opt = sol.get_optimal_weights(max)            % Calculate the optimal weights of the selectable measurements
quality_opt = sol.get_quality(w_opt)            % Calculate quality resulting from optimal measurements
w_subopt = [ones(max, 1); zeros(n-max, 1)]      % Suboptimal weights
quality_subopt = sol.get_quality(w_subopt)      % Calculate quality resulting from suboptimal measurements


%% Estimate model parameter from accomplished measurements
m = 5;                                                                          % Number of accomplished measurements
t_fix = t_opt;                                                                  % Accomplished measurements
v_fix = v_var(w_opt);                                                           % Variances of measurement results at these measurements
eta = model_util.get_fictitious_measurement_results(model, p, t_fix, v_fix);    % Measurement results of the accomplished measurements
sol.set_accomplished_measurements(t_fix, v_fix, eta);                           % Pass accomplished measurements to the solver object
p_lb = [0; 0];                                                                  % Lower bounds of model parameters
p_ub = [10; 5];                                                                 % Upper bounds of model parameters
p_opt = sol.get_optimal_parameters(p_lb, p_ub)                                  % Optimize model parameter from accomplished measurements


%% Calculate gain of additional measurements
sol.set_initial_parameter_estimation(p_opt);   % Update parameter estimation
w_opt = sol.get_optimal_weights(max)';         % Calculate the optimal weights of the selectable measurements
quality_old = sol.get_quality(zeros(n, 1))     % Calculate quality without additional measurements
quality_new = sol.get_quality(w_opt)           % Calculate quality resulting from optimal additional measurements


%% Calculate optimal measurements with constraints
% We are constraining the choice of measurements in such a way that
% distance between two chosen measurements has to be at least 7.5.
A = diag(ones(n, 1)) + diag(ones(n-1, 1), 1) + diag(ones(n-2, 1), 2)    % Matrix for the constraints of the measurements
b = ones(n, 1)                                                          % Vector for the constraints of the measurements
t_opt = sol.get_optimal_measurements(A, b)                              % Calculate the optimal measurements of the selectable measurements considering the constraints

