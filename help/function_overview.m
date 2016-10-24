%% Function Reference
% The classes of the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> are listed here.

%% Solver
% * <matlab:doc('solver') |solver|>           - allows to calculate and optimize the quality of measurements of an experiment and allows to calculate a parameter estimation resulting from accomplished measurements.
% * <matlab:doc('solver_options') |solver_options|>   - represents the options for a solver object.
% * <matlab:doc('solver_edo_options') |solver_edo_options|>                - represents the options for the solver of the experimental design optimization problem.
% * <matlab:doc('solver_po_options') |solver_po_options|>                - represents the options for the solver of the parameter optimization problem.
%

%% Criteria
% * <matlab:doc('criterion') |criterion|>                     - represents an interface for a quality criterion.
% * <matlab:doc('criterion_A') |criterion_A|>                   - is the quality criterion that uses the average variance as the quality.
%

%% Models
% * <matlab:doc('optimal_experimental_design_toolbox.model') |model|>                         - represents an interface for a model.
% * <matlab:doc('model_explicit') |model_explicit|>                - implements the model interface and provides the function value and the first and second derivatives with respect to the parameters of an explicitly given model function.
% * <matlab:doc('model_ivp') |model_ivp|>                     - implements the model interface and provides the solution of an initial value problem and its first and second derivatives with respect to the parameters.
% * <matlab:doc('model_fd') |model_fd|>              - implements the model interface and provides the first and second derivatives with respect to the parameters by finite differences approximations.

%% Utilities
% * <matlab:doc('optimal_experimental_design_toolbox.util') |util|>                          - represents an utility class with static auxiliary methods for the toolbox.
%

%% Application examples
% * <matlab:doc('model_C2') |model_C2|>                      - models the concentration C of the suspended sediment above the marsh surface with two parameters.
% * <matlab:doc('model_C3') |model_C3|>                      - models the concentration C of the suspended sediment above the marsh surface with three parameters.
% * <matlab:doc('model_util') |model_util|>                      - represents an utility class with static auxiliary methods for working with model classes.