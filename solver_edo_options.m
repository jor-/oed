classdef solver_edo_options < handle
% SOLVER_EDO_OPTIONS represents the options for the solver of the experimental design optimization problem.
%
% SOLVER_EDO_OPTIONS Methods:
%     SET_OPTION - changes an option.
%     GET_OPTION - returns the value of an option.
%     GET_SOLVER_ALGORITHM - returns the algorithm to be used to solve the experimental design optimization problem.
%     USE_ALGORITHM_LOCAL_SQP - returns whether to use the local SQP algorithm to solve the relaxed experimental design optimization problem or not.
%     USE_ALGORITHM_DIRECT - returns whether to use the direct algorithm to solve the experimental design optimization problem or not.
%     GET_MAX_FUN_EVALS - returns the maximal model evaluations done by the local SQP solver algorithm.
%     GET_MAX_ITER - returns the maximal iterations done by the local SQP solver algorithm.

%{
---------------------------------------------------------------------------
    Author: Joscha Reimer, jor@informatik.uni-kiel.de, 2010-2013

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
        algorithm_id = 'edo_algorithm';
        algorithm_direct = 'direct';
        algorithm_local_sqp = 'local_sqp';
        
        max_fun_evals_id = 'edo_max_fun_evals';
        max_iter_id = 'edo_max_iter';
        
        error_id__set_option__unknown_option_name = solver_edo_options.get_message_identifier('set_option', 'unknown_option_name');
    end
    
    methods (Access = public)
        
        function this = solver_edo_options(varargin)
        % SOLVER_EDO_OPTIONS creates a SOLVER_EDO_OPTIONS object.
        %
        % Example:
        %     OBJ = SOLVER_EDO_OPTIONS('OPTION1',VALUE1,'OPTION2',VALUE2,...)
        %
        % Input:
        %     'edo_algorithm': the method to be used to solve the experimental
        %         design optimization problem (possible values: 'direct', 
        %         'local_sqp', default: 'local_sqp')
        %     'edo_max_fun_evals': the number of maximal model evaluations done
        %         by the 'local_sqp' solver (possible values: a non-negative integer,
        %         default: 10^3)
        %     'edo_max_iter': the number of maximal iterations of the 'local_sqp' solver
        %         (possible values: a  non-negative integer, default: 10^3)
        %
        % Output:
        %     OBJ: a SOLVER_EDO_OPTIONS object with the passed configurations
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
        % see also SET_OPTION
        %
            
            % set default options
            this.set_option(this.algorithm_id, this.algorithm_local_sqp);
            this.set_option(this.max_fun_evals_id, 10^3);
            this.set_option(this.max_iter_id, 10^3);
            
            % insert passed options
            if mod(nargin, 2) == 0
                for i=1:2:nargin
                    this.set_option(varargin{i}, varargin{i+1});
                end
            else
            	error(this.get_message_identifier('solver_edo_options', 'wrong_number_of_arguments'), 'The number of input arguments is odd. Please check the input arguments.');
            end            
            
        end
        
                
        function set_option(this, name, value)
        % SET_OPTION changes an option.
        %
        % Example:
        %     SOLVER_EDO_OPTIONS_OBJECT.SET_OPTION(NAME, VALUE)
        %
        % Input:
        %     NAME: the name of the option to be changed
        %     VALUE: the new value of the option
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS
        %
            
            % check option name
            if ~ ischar(name)
                error(this.get_message_identifier('set_option', 'name_no_string'), 'The optione name has to be a string.');
            end
            
            % check option value
            switch name
                case this.algorithm_id
                    if ~ (isequal(value, this.algorithm_local_sqp) || ...
                          isequal(value, this.algorithm_direct))
                        error(this.get_message_identifier('set_option', 'unknown_algorithm'), ['The value for "' name '" has to be ' this.algorithm_local_sqp, ' or ' this.algorithm_direct '.']);
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
        %     SOLVER_EDO_OPTIONS_OBJECT.GET_OPTION(NAME)
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
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
            try
                option = this.options.(name);
            catch exception
                error(this.get_message_identifier('get_option', 'unknown_option_name'), ['The option "', name, '" is not supported.']);  
            end
        end
                
        
        function solver_algorithm = get_solver_algorithm(this)
        % GET_SOLVER_ALGORITHM returns the algorithm to be used to solve the experimental design optimization problem.
        %
        % Example:
        %     SOLVER_EDO_OPTIONS_OBJECT.GET_SOLVER_ALGORITHM()
        %
        % Output:
        %     SOLVER_ALGORITHM: the algorithm to be used to solve the
        %                       experimental design optimization problem
        %
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
           solver_algorithm =  this.get_option(this.algorithm_id);
        end
            
        function boolean = use_algorithm_local_sqp(this)
        % USE_ALGORITHM_LOCAL_SQP returns whether to use the local SQP algorithm to solve the relaxed experimental design optimization problem or not.
        %
        % Example:
        %     SOLVER_EDO_OPTIONS_OBJECT.USE_ALGORITHM_LOCAL_SQP()
        %
        % Output:
        %     BOOLEAN: whether to use the local SQP algorithm to solve the
        %         relaxed experimental design optimization problem or not
        %
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
            boolean = isequal(this.get_solver_algorithm(), this.algorithm_local_sqp);
        end
        
        function boolean = use_algorithm_direct(this)
        % USE_ALGORITHM_DIRECT returns whether to use the direct algorithm to solve the experimental design optimization problem or not.
        %
        % Example:
        %     SOLVER_EDO_OPTIONS_OBJECT.USE_ALGORITHM_DIRECT()
        %
        % Output:
        %     BOOLEAN: whether to use the direct algorithm to solve the
        %         experimental design optimization problem or not
        %
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
            boolean = isequal(this.get_solver_algorithm(), this.use_algorithm_direct);
        end
        
        function max_fun_evals = get_max_fun_evals(this)
        % GET_MAX_FUN_EVALS returns the maximal model evaluations done by the local SQP solver algorithm.
        %
        % Example:
        %     SOLVER_EDO_OPTIONS_OBJECT.GET_MAX_FUN_EVALS()
        %
        % Output:
        %     MAX_FUN_EVALS: the maximal model evaluations done by the 
        %                    local SQP solver algorithm
        %
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
            max_fun_evals = this.get_option(this.max_fun_evals_id);
        end
        
        function max_iter = get_max_iter(this)
        % GET_MAX_ITER returns the maximal iterations done by the local SQP solver algorithm.
        %
        % Example:
        %     SOLVER_EDO_OPTIONS_OBJECT.GET_MAX_ITER()
        %
        % Output:
        %     MAX_ITER: the maximal iterations done by the 
        %               local SQP solver algorithm
        %
        % see also SOLVER_EDO_OPTIONS.SOLVER_EDO_OPTIONS, SET_OPTION
        %
        
            max_iter = this.get_option(this.max_iter_id);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function s = get_message_identifier(method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in methods of these object.
        %
        % Example:
        %     ID = SOLVER_EDO_OPTIONS.GET_MESSAGE_IDENTIFIER(METHOD, MNEMONIC)
        %
        % Input:
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
        
            s = util.get_message_identifier('solver_edo_options', method, mnemonic);
        end
        
    end
end

