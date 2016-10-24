%% User's guide
% This user's guide demonstrates the features of the
% <matlab:doc('optimal_experimental_design_toolbox') |Optimal Experimental Design Toolbox|>
% and their application.



%% Specify the model
% The first step to optimizing of model parameters or corresponding
% measurement conditions with the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> is to specify the model function. This is
% done by providing an object that implements the <matlab:doc('model')
% |model|> interface, which prescribes methods for the model function and
% its first and second derivative with respect to the model parameters.

model_object = model_class();

%%
% If the model function is given in an explicit form you can use the
% <matlab:doc('model_explicit') |model_explicit|> class. The necessary
% function values and derivatives will then be calculated automatically.

model_object = model_explicit(model_function, parameter_variables, experimental_design_variables);
%%
% If the model function is given as an initial value problem you can use
% the <matlab:doc('model_ivp') |model_ivp|> class. There only the initial
% value problem must be specified. The necessary function values and
% derivatives will then be calculated automatically.

model_object = model_ivp(differential_equation, parameter_variables, ...
    independent_variable, initial_value, dependent_variable, integration_interval);

%%
% If you are using the <matlab:doc('model_ivp') |model_ivp|> class you can
% obtain a speedup by calculating the necessary derivations of the model
% in parallel. To achieve that you only have to open a |matlabpool| before
% the calulations. 

matlabpool open number_of_cores

%%
% If you have an arbitrary implementation of the model function without the
% derivatives with respect to the parameters you can implement the
% <matlab:doc('model_fd') |model_fd|> interface instead of the
% <matlab:doc('model') |model|> interface. This provides the first and
% second derivatives with respect to the parameters by finite differences
% approximations.






%% Create and configure an |solver| object
% The next step is to create and configure an <matlab:doc('solver')
% |solver|> object. First of all you have to create an 
% <matlab:doc('solver') |solver|> object via the
% <matlab:doc('solver.solver') |constructor|>.

solver_object = solver(model_object, ... 
    parameter_estimation, measurement_points, corresponding_variances);

%%
% Then you have to pass the previously created <matlab:doc('model')
% |model|> object presenting your model function to your
% <matlab:doc('solver') |solver|> object. This can be done via the
% <matlab:doc('solver.solver') |constructor|> or
% the <matlab:doc('solver.set_model') |set_model|> method.

solver_object.set_model(model_object);

%%
% Further an approximation of the parameter values must be specified via
% the <matlab:doc('solver.solver') |constructor|>
% or the <matlab:doc('solver.set_initial_parameter_estimation')
% |set_initial_parameter_estimation|> method. This can arise for example
% from experience or previous measurements.

solver_object.set_initial_parameter_estimation(parameter_estimation);

%%
% Next you have to specify selectable measurements and there
% corresponding variances of the measurement errors. This can be done via
% the <matlab:doc('solver.solver') |constructor|>
% and the <matlab:doc('solver.set_selectable_measurements')
% |set_selectable_measurements|> method.

solver_object.set_selectable_measurements(measurement_points, ...
    corresponding_variances);

%%
% If you already have done previous measurements, this can be included. The
% measurements and there corresponding variances of the measurement
% errors can be set via the
% <matlab:doc('solver.set_accomplished_measurements')
% |set_accomplished_measurements|> method.

solver_object.set_accomplished_measurements(measurement_points, ...
    corresponding_variances, measurement_results);





%% Optimize measurements
% Now you can calculate the optimal measurements by using the
% <matlab:doc('solver.get_optimized_weights')
% |get_optimized_weights|> method of your solver object.
% The result of the <matlab:doc('solver.get_optimized_weights')
% |get_optimized_weights|> method is a weight 0 or 1 for each
% measurements. 1 means the measurement has been chosen and 0 means
% the measurement has not been chosen. You can pass linear constraints
% to the measurements. Thereby, for example, a restriction to the cost
% of the measurements can be achieved.

weights = solver_object.get_optimized_weights(A, b);

%%
% A common case is that the number of selected measurements should be
% restricted. This can be achieved directly as follows. 

weights = solver_object.get_optimized_weights(maximal);





%% Quality of measurements
% You can also calculate the quality of a set of measurements with the
% <matlab:doc('solver.get_quality') |get_quality|> method. This allows, for
% example, to calculate the quality resulting by omitting some measurements.
    
quality = solver_object.get_quality(weights);





%% Optimize model parameters
% If you have done previous measurements you can optimize the model
% parameters with your measurement results. For that you can use the
% <matlab:doc('solver.get_optimal_parameters') |get_optimal_parameters|>
% method.

parameter = solver_object.get_optimal_parameters(lower_bound, upper_bound);





%% Change options
% Several options are available to configure the optimization of the
% experimental designs and model parameters according to your needs. You
% can set your desired option with the <matlab:doc('solver.set_option')
% |set_option|> method of your <matlab:doc('solver') |solver|> object.

solver_object.set_option('option1', value1);
solver_object.set_option('option2', value2);

%%
% The possible options are
%
% * |'parameter_estimation'|: whether a parameter estimation should be
% performed befor the optimal design estimation or not (_possible values:_
% |'yes'|, |'no'|, _default:_ |'no'|)
%
% * |'estimation_method'|: the method of the estimation of the quality of
% the experimental design (_possible values:_ |'point'|, |'region'|,
% _default:_ |'region'|)
%
% * |'alpha'|: the confidence level for the region estimation (_possible
% value:_ a scalar with |0 < alpha < 1|, _default:_ |0.95|) 
%
% * |'debug'|: the level of debug information to be displayed (_possible
% values:_ a non-negative integer, _default:_ |0| (no debug informations))
%
% * |'criterion'|: the criterion for the quality of an experimental design
% (_possible values:_ an object of the <matlab:doc('criterion') |criterion|>
% class, _default:_ a <matlab:doc('criterion_A') |criterion_A|> object)
%
% * |'edo_algorithm'|: the method to be used to solve the experimental design
% optimization problem (_possible values_: |'direct'|, |'local_sqp'|,
% _default_: |'local_sqp'|)
%
% * |'edo_scale_parameters'|:  whether the parameters should be scaled
% or not (_possible values:_ |'yes'|, |'no'|, _default:_ |'yes'|)
%
% * |'edo_max_fun_evals'|: the number of maximal model evaluations done by the
% 'local_sqp' solver (_possible values_: a non-negative integer, _default_: |10^3|)
%
% * |'edo_max_iter'|: the number of maximal iterations of the 'local_sqp' solver
% (_possible values_: a  non-negative integer, _default_: |5*10^1|) 
%
% * |'po_algorithm'|: the method to be used to solve the parameter
% optimization problem (_possible values_: |'trust-region-reflective'|, 
% |'levenberg-marquardt'|, _default_: |'trust-region-reflectiv'|)
%
% * |'po_scale_parameters'|: whether the parameters have to be scaled
% for the optimization (_possible values_: |'yes'|, |'no'|,
%  _default_: |'yes'|)
%
% * |'po_max_fun_evals'|: the number of maximal model evaluations done
% by the solver (_possible values_: a non-negative integer, 
% _default_: |3 * 10^3|)
%
% * |'po_max_iter'|: the number of maximal iterations of the solver
% (_possible values_: a  non-negative integer, _default_: |3 * 10^3|)
%
