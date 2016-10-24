%% Product overview 
% The <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> allows to optimize model parameters and
% corresponding measurement conditions.


%% Introduction
% Often models contain roughly known model parameters, which should to be
% determined more accurately. These parameters can be optimized, so that
% the model is as consistent as possible with results of measurements.
% These measurement results are obtained under different measurement
% conditions, which crucially affect the information content of the
% measurement results. These measurement conditions can also be optimized
% so that the information contain of the corresponding measurement results is
% maximized. As a result fewer measurements have to be carried out to gain
% sufficient accurate model parameters which saves time, money and effort.


%% Optimize measurement conditions
% You only have to specify the model function, the selectable measurements
% and the variances of the corresponding measurement errors. Than the
% <matlab:doc('optimal_experimental_design_toolbox') |Optimal Experimental Design Toolbox|>
% can calculate which measurements you have to use to be able to
% determine the parameters optimally. Optimally here means that the
% expected error, which will be done by determining the parameters, is
% minimal.


%% Constrain the measurements
% In addition it is possible to define linear constraints to the measurement
% points. Thus, for example, the total cost of the measurement that will be
% chosen can be limited. A simple constraint, for example, can be a maximal
% number of measurements.


%% Take previous measurements into account 
% If for the modeled value measurement results already exist, the corresponding
% measurements and variances can also be considered. The |Optimal Experimental
% Design Toolbox| can tell you how great the benefit of additional
% measurements would be and thereby, if further measurements are necessary.


%% Optimize model parameters
% Furthermore the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> can optimize the parameters of models
% by means of accomplished measurements. Box constraints on the model
% parameters can be considered in the optimization. 
