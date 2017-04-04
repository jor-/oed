classdef util
    
% UTIL represents an utility class with static auxiliary methods for the toolbox.
%
% UTIL Methods:
%   MAKE_COLUMN_VECTOR - makes the passed vector to a column vector.
%   MAKE_ROW_VECTOR - makes the passed vector to a row vector.
%   MAKE_SYM - returns the parsed string as symbolic formula.
%   GET_MESSAGE_IDENTIFIER - returns the identifier for an error or a warning raised in this toolbox.
        
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

    properties
    end
    
    methods (Access = public, Static)
        
        function v = make_column_vector(v)
        % MAKE_COLUMN_VECTOR makes the passed vector to a column vector.
        %
        % Example:
        %     V = UTIL.MAKE_COLUMN_VECTOR(V)
        %
        % Input:
        %     V: the vector
        %
        % Output:
        %     V: the vector as a column vector
        %
        % Throws:
        %     An error if V is not a vector.
        %
        
            sizeV = size(v);
            if sizeV(2) > 1
                if sizeV(1) > 1
                    display(v);
                    error(util.get_message_identifier('util', 'make_row_vector', 'vector_expected'), 'vector expected, not a matrix');
                else
                    v = v.';
                end
            end
        end
        
        function v = make_row_vector(v)
        % MAKE_ROW_VECTOR makes the passed vector to a row vector.
        %
        % Example:
        %     V = UTIL.MAKE_ROW_VECTOR(V)
        %
        % Input:
        %     V: the vector
        %
        % Output:
        %     V: the vector as a row vector
        %
        % Throws:
        %     An error if V is not a vector.
        %
        
            sizeV = size(v);
            if sizeV(1) > 1
                if sizeV(2) > 1
                    display(v);
                    error(util.get_message_identifier('util', 'make_row_vector', 'vector_expected'), 'vector expected, not a matrix');
                else
                    v = v.';
                end
            end
        end
        
        function f_sym = make_sym(f)
        % MAKE_SYM returns the parsed string or cell array of strings as symbolic formulas.
        %
        % Example:
        %     F_SYM = UTIL.MAKE_SYM(F)
        %
        % Input:
        %     F: the string or the cell array of strings
        %
        % Output:
        %     F_SYM: the symbolic representation of F
        %
        % Throws:
        %     An error if F is not a string, a cell array of strings or a
        %     symbolic function.
        %
        
            if ischar(f) || iscellstr(f)
                f_sym = sym(f);
            elseif isa(f, 'sym')
                f_sym = f;
            else
                error(util.get_message_identifier('util', 'make_sym', 'wrong__input_type'), ['The input should be a string or a symbolic function. But it is a ' class(f)]);
            end
        end
        
        function Jacobian = approximate_Jacobian(F, x, order, Fx)
            if nargin < 3
                order = 2;
            end
            
            n = length(x);
            h = util.get_h_for_finite_differences(x);
            
            switch order
                case 1
                    if nargin < 4
                        Fx = F(x);
                    end
                    for i = 1:n
                        e = zeros(size(x));
                        e(i) = 1;

                        Jacobian(:, i) = util.make_column_vector((F(x + h(i) * e) - Fx) / (2*h(i)));
                    end
                case 2
                    for i = 1:n
                        e = zeros(size(x));
                        e(i) = 1;

                        Jacobian(:, i) = util.make_column_vector((F(x + h(i) * e) - F(x - h(i) * e)) / (2*h(i)));
                    end
                otherwise
                    error(util.get_message_identifier('util', 'approximate_Jacobian', 'order not implemented'), 'Methods are only implemented for order 1 and 2.');
            end
        end
        
        
        function Hessian = approximate_Hessian(F, x, grad_Fx)
            grad_F = @(x) (util.approximate_Jacobian(F, x, 2));
            if nargin < 3
                grad_Fx = grad_F(x);
            end
            
            Hessian = util.approximate_Jacobian(grad_F, x, 1, grad_Fx);
            
            Hessian = (Hessian + Hessian') / 2;
        end
        
        function id = get_message_identifier(class, method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in this toolbox.
        %
        % Example:
        %     ID = UTIL.GET_MESSAGE_IDENTIFIER(CLASS, METHOD, MNEMORIC)
        %
        % Input:
        %     CLASS: the class in which an error or a warning occurred
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
        
            id = ['optimal_experimental_design_toolbox:', class, ':', method, ':', mnemonic];
        end
                                
    end  
    
    methods (Access = protected, Static)
        function h = get_h_for_finite_differences(x)
            typ = sqrt(eps);
            n = length(x);           
            h = zeros(n, 1);

            for i = 1:n               
                if x(i) == 0
                    h(i) = typ;
                else
                    h(i) = max(abs(x(i)), typ) * sign(x(i));
                end
            end
            h = h * sqrt(eps) * 10^2;
        end
    end
end

