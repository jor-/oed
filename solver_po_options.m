classdef solver_po_options < handle
% SOLVER_PO_OPTIONS represents the options for the solver of the parameter optimization problem.
%
% SOLVER_PO_OPTIONS Methods:
%     SOLVER_PO_OPTIONS - creates a SOLVER_PO_OPTIONS object.
%     SET_OPTION - changes an option.
%     GET_OPTION - returns the value of an option.
%     GET_SOLVER_ALGORITHM - returns the algorithm to be used to solve the parameter optimization problem.
%     USE_ALGORITHM_TRUST_REGION_REFLECTIVE - returns whether to use the Trust-Region-Reflective algorithm to solve the parameter optimization problem or not.
%     USE_ALGORITHM_LEVENBERG_MARQUARDT - returns whether to use the Levenberg-Marquardt algorithm to solve the parameter optimization problem or not.
%	  SCALE_PARAMETERS - returns whether the parameters have to be scaled for the optimization or not.
%     GET_MAX_FUN_EVALS - returns the maximal model evaluations done by the solver algorithm.
%     GET_MAX_ITER - returns the maximal iterations done by the solver algorithm.
%     GET_MESSAGE_IDENTIFIER - returns the identifier for an error or a warning raised in methods of these object.

%{
---------------------------------------------------------------------------
    Copyright (C) 2010-2013 Joscha Reimer jor@informatik.uni-kiel.de

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
    end
    
    properties (Constant)        
        algorithm_id = 'po_algorithm';
        algorithm_trust_region_reflective = 'trust-region-reflective';
        algorithm_levenberg_marquardt = 'levenberg-marquardt';
        
        scale_parameter_id = 'po_scale_parameter';
        scale_parameter_yes = 'yes';
        scale_parameter_no = 'no';
        
        max_fun_evals_id = 'po_max_fun_evals';
        max_iter_id = 'po_max_iter';
        
        error_id__set_option__unknown_option_name = solver_po_options.get_message_identifier('set_option', 'unknown_option_name');
    end
    
    methods (Access = public)
        
        function this = solver_po_options(varargin)
        % SOLVER_PO_OPTIONS creates a SOLVER_PO_OPTIONS object.
        %
        % Example:
        %     OBJ = SOLVER_PO_OPTIONS('OPTION1',VALUE1,'OPTION2',VALUE2,...)
        %
        % Input:
        %     'po_algorithm': the method to be used to solve the parameter
        %         optimization problem (possible values: 'trust-region-reflective', 
        %         'levenberg-marquardt', default: 'trust-region-reflectiv')
        %     'po_scale_parameter': whether the parameter have to be scaled
        %         for the optimization (possible values: 'yes', 'no',
        %         default: 'yes')
        %     'po_max_fun_evals': the number of maximal model evaluations done
        %         by the solver (possible values: a non-negative integer,
        %         default: 3 * 10^3)
        %     'po_max_iter': the number of maximal iterations of the solver
        %         (possible values: a  non-negative integer, default: 3 * 10^3)
        %
        % Output:
        %     OBJ: a SOLVER_PO_OPTIONS object with the passed configurations
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
        % see also SET_OPTION
        %
            
            % set default options
            this.set_option(this.algorithm_id, this.algorithm_trust_region_reflective);
            this.set_option(this.scale_parameter_id, this.scale_parameter_yes);
            this.set_option(this.max_fun_evals_id, 3 * 10^3);
            this.set_option(this.max_iter_id, 3 * 10^3);
            
            % insert passed options
            if mod(nargin, 2) == 0
                for i=1:2:nargin
                    this.set_option(varargin{i}, varargin{i+1});
                end
            else
            	error(this.get_message_identifier('solver_po_options', 'wrong_number_of_arguments'), 'The number of input arguments is odd. Please check the input arguments.');
            end            
            
        end        
        
        
        function set_option(this, name, value)
        % SET_OPTION changes an option.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.SET_OPTION(NAME, VALUE)
        %
        % Input:
        %     NAME: the name of the option to be changed
        %     VALUE: the new value of the option
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS
        %
            
            % check option name
            if ~ ischar(name)
                error(this.get_message_identifier('set_option', 'name_no_string'), 'The optione name has to be a string.');
            end
            
            % check option value
            switch name
                case this.algorithm_id
                    if ~ (isequal(value, this.algorithm_trust_region_reflective) || ...
                          isequal(value, this.algorithm_levenberg_marquardt))
                        error(this.get_message_identifier('set_option', 'unknown_algorithm'), ['The value for "' name '" has to be ' this.algorithm_local_sqp, ' or ' this.algorithm_direct '.']);
                    end
                case this.scale_parameter_id
                    if ~ (isequal(value, this.scale_parameter_yes) || isequal(value, this.scale_parameter_no))
                        error(this.get_message_identifier('set_option', 'unknown_scale_parameter_option'),  ['The value for "', name, '" has to be ', this.scale_parameter_yes, ' or ', this.scale_parameter_no, '.']);
                    end
                case this.max_fun_evals_id
                    if ~ (isscalar(value) && value == fix(value))
                        error(this.get_message_identifier('set_option', 'no_scalar_integer'), ['The value for "', name, '" has to be scalar integer.']);
                    end
                case this.max_iter_id
                    if ~ (isscalar(value) && value == fix(value))
                        error(this.get_message_identifier('set_option', 'no_scalar_integer'), ['The value for "', name, '" has to be scalar integer.']);
                    end
                otherwise
                    error(this.error_id__set_option__unknown_option_name, ['The option "', name, '" is not supported.']); 
            end
            
            % update option
            this.options.(name) = value;            
        end        
        
        
        function option = get_option(this, name)
        % GET_OPTION returns the value of an option.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.GET_OPTION(NAME)
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
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
            try
                option = this.options.(name);
            catch exception
                error(this.get_message_identifier('get_option', 'unknown_option_name'), ['The option "', name, '" is not supported.']);  
            end
        end        
        
        
        function solver_algorithm = get_solver_algorithm(this)
        % GET_SOLVER_ALGORITHM returns the algorithm to be used to solve the parameter optimization problem.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.GET_SOLVER_ALGORITHM()
        %
        % Output:
        %     SOLVER_ALGORITHM: the algorithm to be used to solve the
        %                       paramater optimization problem
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
           solver_algorithm =  this.get_option(this.algorithm_id);
        end
            
        function boolean = use_algorithm_trust_region_reflective(this)
        % USE_ALGORITHM_TRUST_REGION_REFLECTIVE returns whether to use the Trust-Region-Reflective algorithm to solve the parameter optimization problem or not.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.USE_ALGORITHM_TRUST_REGION_REFLECTIVE()
        %
        % Output:
        %     BOOLEAN: whether to use the local SQP algorithm to solve the
        %         parameter optimization problem or not
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
            boolean = isequal(this.get_solver_algorithm(), this.algorithm_trust_region_reflective);
        end
        
        function boolean = use_algorithm_levenberg_marquardt(this)
        % USE_ALGORITHM_LEVENBERG_MARQUARDT returns whether to use the Levenberg-Marquardt algorithm to solve the parameter optimization problem or not.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.USE_ALGORITHM_LEVENBERG_MARQUARDT()
        %
        % Output:
        %     BOOLEAN: whether to use the Levenberg-Marquardt algorithm to solve the
        %         parameter optimization problem or not
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
            boolean = isequal(this.get_solver_algorithm(), this.algorithm_levenberg_marquardt);
        end
        
        function boolean = scale_parameters(this)
        % SCALE_PARAMETERS returns whether the parameters have to be scaled for the optimization or not.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.USE_ALGORITHM_LEVENBERG_MARQUARDT_SCALED()
        %
        % Output:
        %     BOOLEAN: whether the parameters have to be scaled for the
        %              optimization or not.
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
            boolean = isequal(this.get_option(this.scale_parameter_id), this.scale_parameter_yes);
        end
        
        function max_fun_evals = get_max_fun_evals(this)
        % GET_MAX_FUN_EVALS returns the maximal model evaluations done by the solver algorithm.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.GET_MAX_FUN_EVALS()
        %
        % Output:
        %     MAX_FUN_EVALS: the maximal model evaluations done by the 
        %                    solver algorithm
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
            max_fun_evals = this.get_option(this.max_fun_evals_id);
        end
        
        function max_iter = get_max_iter(this)
        % GET_MAX_ITER returns the maximal iterations done by the solver algorithm.
        %
        % Example:
        %     SOLVER_PO_OPTIONS_OBJECT.GET_MAX_ITER()
        %
        % Output:
        %     MAX_ITER: the maximal iterations done by the 
        %               solver algorithm
        %
        % see also SOLVER_PO_OPTIONS.SOLVER_PO_OPTIONS, SET_OPTION
        %
        
            max_iter = this.get_option(this.max_iter_id);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function s = get_message_identifier(method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in methods of these object.
        %
        % Example:
        %     ID = SOLVER_PO_OPTIONS.GET_MESSAGE_IDENTIFIER(METHOD, MNEMONIC)
        %
        % Input:
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
        
            s = util.get_message_identifier('solver_po_options', method, mnemonic);
        end
        
    end
end

