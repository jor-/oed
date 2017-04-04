classdef solver_options < handle
% SOLVER_OPTIONS represents the options for a solver object.
%
% SOLVER_OPTIONS Methods:
%   SET_OPTION - changes an option.
%   GET_OPTION - returns the value of an option.
%   GET_ALPHA - returns the confidence level for the region estimation.
%   GET_CRITERION - returns the criterion for the quality of an experimental design.
%   GET_ESTIMATION_METHOD - returns the method to use to estimate the quality of an experimental design.
%   USE_DEBUG - returns whether debug informations will be output or not.
%   USE_ESTIMATION_METHOD_REGION - returns whether to use the region method to estimate the quality of an experimental design or not.
%   USE_PARAMETER_ESTIMATION - returns whether to perform a parameter estimation before the optimization of the experimental design or not.
%   GET_MESSAGE_IDENTIFIER - returns the identifier for an error or a warning raised in methods of these object.
%   GET_SOLVER_EDO_OPTIONS - returns the options for the solver of the experimental design optimization problem.
%   GET_SOLVER_PO_OPTIONS - returns the options for the solver of the parameter optimization problem.

%{
---------------------------------------------------------------------------
    Copyright (C) 2010-2017 Joscha Reimer jor@informatik.uni-kiel.de

    This file is part of the Optimal Experimental Design Toolbox.

    The Optimal Experimental Design Toolbox is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    The Optimal Experimental Design Toolbox is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Optimal Experimental Design Toolbox. If not, see
    <http://www.gnu.org/licenses/>.
---------------------------------------------------------------------------
%}
    
    properties (Access = private)
        options = struct;
        
        criterion_class_name = 'criterion';
        solver_edo_options_class_name = 'solver_edo_options';
        solver_po_options_class_name = 'solver_po_options';
    end
        
    properties (Constant)
        debug_id = 'debug';
        alpha_id = 'alpha';
        
        criterion_id = 'criterion';
        
        estimation_method_id = 'estimation_method';
        estimation_method_point = 'point';
        estimation_method_region = 'region';       
        
        parameter_estimation_id = 'parameter_estimation';
        parameter_estimation_yes = 'yes';
        parameter_estimation_no = 'no';
        
        solver_edo_options_id = 'solver_edo_options';
        solver_po_options_id = 'solver_po_options';
    end
    
    methods (Access = public)
        
        function this = solver_options(varargin)
        % SOLVER_OPTIONS creates a SOLVER_OPTIONS object.
        %
        % Example:
        %     OBJ = SOLVER_OPTIONS('OPTION1',VALUE1,'OPTION2',VALUE2,...)
        %
        % Input:
        %     'parameter_estimation': whether a parameter estimation should
        %         be performed before the optimal design estimation or not
        %         (possible values: 'yes', 'no', default: 'no')
        %     'estimation_method': the method of the estimation of the
        %         quality of the experimental design (possible values: 
        %         'point', 'region', default: 'region')
        %     'alpha': the confidence level for the region estimation 
        %         (possible value: a scalar with 0 < alpha < 1,
        %         default: 0.95)
        %     'debug': the level of debug information to be displayed
        %         (possible value: a non-negative integer,
        %         default: 0 (no debug informations))
        %     'criterion': the criterion for the quality of an experimental
        %         design (possible values: an object of the CRITERION class,
        %         default: a CRITERION_A object)
        %     'solver_edo_options': the solver to be used to solve the
        %         experimental design optimization problem (possible
        %         values: an object of the SOLVER_EDO_OPTIONS class,
        %         default: a default SOLVER_EDO_OPTIONS object)
        %     'solver_po_options': the solver to be used to solve the
        %         parameter optimization problem (possible
        %         values: an object of the SOLVER_PO_OPTIONS class,
        %         default: a default SOLVER_PO_OPTIONS object)
        %
        % Output:
        %     OBJ: a SOLVER_OPTIONS object with the passed configurations
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
        % see also CRITERION, CRITERION_A, SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
            % set default options
            this.set_option(this.debug_id, 0);
            this.set_option(this.alpha_id, 0.95);
            this.set_option(this.criterion_id, criterion_A());
            this.set_option(this.estimation_method_id, this.estimation_method_region);
            this.set_option(this.parameter_estimation_id, this.parameter_estimation_no);
            this.set_option(this.solver_edo_options_id, solver_edo_options());
            this.set_option(this.solver_po_options_id, solver_po_options());
            
            % insert passed options
            if mod(nargin, 2) == 0
                for i=1:2:nargin
                    this.set_option(varargin{i}, varargin{i+1});
                end
            else
                error(this.get_message_identifier('solver_options', 'wrong_number_of_arguments'), 'The number of input arguments is odd. Please check the input arguments.');
            end            
            
        end
                
        
        function set_option(this, name, value)
        % SET_OPTION changes an option.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.SET_OPTION(NAME, VALUE)
        %
        % Input:
        %     NAME: the name of the option to be changed
        %     VALUE: the new value of the option
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS, CRITERION, CRITERION_A, SOLVER_EDO_OPTIONS
        %
            
            % check option name
            if ~ ischar(name)
                error(this.get_message_identifier('set_option', 'name_no_string'), 'The optione name has to be a string.');
            end
            
            % check option value
            switch name
                case this.debug_id
                    if ~ isscalar(value)
                        error(this.get_message_identifier('set_option', 'no_scalar'), ['The value for "', name, '" has to be scalar.']);
                    end
                    this.options.(name) = value;
                case this.alpha_id
                    if ~ isscalar(value) || value <= 0 || value >= 1
                        error(this.get_message_identifier('set_option', 'not_between_0_and_1'), ['The value for "', name, '" has to be between 0 and 1.']);
                    end
                    this.options.(name) = value;
                case this.criterion_id
                    if ~ isa(value, this.criterion_class_name)
                        error(this.get_message_identifier('set_option', 'no_criterion_class'), ['The value for "', name, '" has to be a an object of the ', this.criterion_class_name, ' class.']);
                    end
                    this.options.(name) = value;
                case this.estimation_method_id
                    if ~ (isequal(value, this.estimation_method_point) || isequal(value, this.estimation_method_region))
                        error(this.get_message_identifier('set_option', 'unknown_estimation_method_option'), ['The value for "', name, '" has to be ', this.estimation_method_point, ' or ', this.estimation_method_region, '.']);
                    end
                    this.options.(name) = value;
                case this.parameter_estimation_id
                    if ~ (isequal(value, this.parameter_estimation_yes) || isequal(value, this.parameter_estimation_no))
                        error(this.get_message_identifier('set_option', 'unknown_parameter_estimation_option'),  ['The value for "', name, '" has to be ', this.parameter_estimation_yes, ' or ', this.parameter_estimation_no, '.']);
                    end
                    this.options.(name) = value;
                case this.solver_edo_options_id
                    if ~ isa(value, this.solver_edo_options_class_name)
                        error(this.get_message_identifier('set_option', 'no_solver_edo_option_class'), ['The value for "', name, '" has to be a subclass of the ', this.solver_edo_options_class_name, ' class.']);
                    end
                    this.options.(name) = value;
                case this.solver_po_options_id
                    if ~ isa(value, this.solver_po_options_class_name)
                        error(this.get_message_identifier('set_option', 'no_solver_po_option_class'), ['The value for "', name, '" has to be a subclass of the ', this.solver_po_options_class_name, ' class.']);
                    end
                    this.options.(name) = value;
                otherwise
                    try
                        solver_edo_options = this.get_solver_edo_options();
                        solver_edo_options.set_option(name, value);
                    catch exception
                        if isequal(exception.identifier, solver_edo_options.error_id__set_option__unknown_option_name)
                            try
                                solver_po_options = this.get_solver_po_options();
                                solver_po_options.set_option(name, value);
                            catch exception
                                if isequal(exception.identifier, solver_po_options.error_id__set_option__unknown_option_name)
                                    error(this.get_message_identifier('set_option', 'unknown_option_name'), ['The option "', name, '" is not supported.']);
                                else
                                    rethrow(exception);
                                end
                            end
                        else
                            rethrow(exception);
                        end
                    end
            end            
        end        
        
        
        function option = get_option(this, name)
        % GET_OPTION returns the value of an option.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.GET_OPTION(NAME)
        %
        % Input:
        %     NAME: the name of the option which value will be returned
        %
        % Output:
        %     VALUE: the value of the passed option
        %
        % Throws:
        %     An error if the option doesn't exist.
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS, SET_OPTION
        %
        
            try
                option = this.options.(name);
            catch exception
                error(this.get_message_identifier('get_option', 'unknown_option_name'), ['The option "', name, '" is not supported.']);  
            end
        end
        
        
        function alpha = get_alpha(this)
        % GET_ALPHA returns the confidence level for the region estimation.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.GET_ALPHA()
        %
        % Output:
        %     APLHA: the confidence level for the region estimation
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
           alpha =  this.get_option(this.alpha_id);
        end
        
        function criterion = get_criterion(this)
        % GET_CRITERION returns the criterion for the quality of an experimental design.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.GET_CRITERION()
        %
        % Output:
        %     CRITERION: the criterion for the quality of an experimental design
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS, CRITERION
        %
        
           criterion =  this.get_option(this.criterion_id);
        end
        
        function estimation_method = get_estimation_method(this)
        % GET_ESTIMATION_METHOD returns the method to use to estimate the quality of an experimental design.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.GET_ESTIMATION_METHOD()
        %
        % Output:
        %     ESTIMATION_METHOD: the method to use to estimate the quality of an experimental design
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
           estimation_method =  this.get_option(this.estimation_method_id);
        end
        
        function boolean = use_debug(this, debug)
        % USE_DEBUG returns whether debug informations will be output or not.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.USE_DEBUG(DEBUG)
        %
        % Input:
        %     DEBUG: the debug-level which will be checked
        %            (optional, default: 1)
        %
        % Output:
        %     USE_DEBUG: whether the current debug-level is greater or
        %                equal the passed debug-level
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
            if nargin <= 1
                debug = 1;
            end
            boolean = this.get_option(this.debug_id) >= debug;
        end
            
        function boolean = use_estimation_method_region(this)
        % USE_ESTIMATION_METHOD_REGION returns whether to use the region method to estimate the quality of an experimental design or not.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.USE_ESTIMATION_METHOD_REGION()
        %
        % Output:
        %     BOOLEAN: whether to use the region method to estimate the
        %              quality of an experimental design or not
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
            boolean = isequal(this.get_estimation_method(), this.estimation_method_region);
        end
        
        function boolean = use_parameter_estimation(this)
        % USE_PARAMETER_ESTIMATION returns whether to perform a parameter estimation before the optimization of the experimental design or not.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.USE_PARAMETER_ESTIMATION()
        %
        % Output:
        %     BOOLEAN: whether to perform a parameter estimation before the
        %              optimization of the experimental design or not
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
            boolean = isequal(this.get_option(this.parameter_estimation_id), this.parameter_estimation_yes);
        end
        
        function solver_edo_options = get_solver_edo_options(this)
        % GET_SOLVER_EDO_OPTIONS returns the options for the solver of the experimental design optimization problem.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.GET_SOLVER_EDO_OPTIONS()
        %
        % Output:
        %     SOLVER_EDO_OPTIONS: the options for the solver of the experimental design optimization problem
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
           solver_edo_options =  this.get_option(this.solver_edo_options_id);
        end
        
        function solver_po_options = get_solver_po_options(this)
        % GET_SOLVER_PO_OPTIONS returns the options for the solver of the parameter optimization problem.
        %
        % Example:
        %     SOLVER_OPTIONS_OBJECT.GET_SOLVER_PO_OPTIONS()
        %
        % Output:
        %     SOLVER_PO_OPTIONS: the options for the solver of the parameter optimization problem
        %
        % see also SOLVER_OPTIONS.SOLVER_OPTIONS
        %
        
           solver_po_options =  this.get_option(this.solver_po_options_id);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function s = get_message_identifier(method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in methods of these object.
        %
        % Example:
        %     ID = SOLVER_OPTIONS.GET_MESSAGE_IDENTIFIER(METHOD, MNEMONIC)
        %
        % Input:
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
        
            s = util.get_message_identifier('solver_options', method, mnemonic);
        end
        
    end
end

