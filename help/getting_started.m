%% Getting Started
% The requirements of the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> and ways to get further help to the
% <matlab:doc('optimal_experimental_design_toolbox') |Optimal Experimental Design Toolbox|>
% are listed here.


%% System Requirements
% To use some of the features of the <matlab:doc('optimal_experimental_design_toolbox') |Optimal Experimental Design Toolbox|> some
% additional |Matlab Toolboxes| are required.
%
% * |Optimization Toolbox| - The |Optimization Toolbox| is needed for the 
% local SQP solver and the parameter estimation.
%
% * |Statistics Toolbox| - The |Statistics Toolbox| is needed for the region
% estimation.
%
% * |Symbolic Math Toolbox| - The |Symbolic Math Toolbox| is needed for the
% explicit model (<matlab:doc('model_explicit') |model_explicit|>) and the
% initial value problem model (<matlab:doc('model_ivp') |model_ivp|>).
%
% * |Parallel Computing Toolbox| - If the |Parallel Computing Toolbox| is
% available you can achieve a speedup for the initial value problem model
% (<matlab:doc('model_ivp') |model_ivp|>).


%% Tested environment
% The <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|> was tested with the following |Matlab|
% and |Toolbox| versions.
V = ver('Matlab');
Matlab_Version = V.Version
V = ver('optim');
Optimization_Toolbox_Version = V.Version
V = ver('stats');
Statistics_Toolbox_Version = V.Version
V = ver('symbolic');
Symbolic_Math_Toolbox_Version = V.Version
V = ver('distcomp');
Parallel_Computing_Toolbox_Version = V.Version

    

%% Further Help
% For detailed help for individual classes or methods of the
% <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|>, please use the |Matlab|
% <matlab:doc('help') |help|> and <matlab:doc('doc') |doc|> functions.
% Another way is to look directly in the commented source code. Further you
% can use the search function of the |Matlab| help browser to search through
% this help documentation.