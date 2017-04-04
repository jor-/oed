classdef model_explicit < model
% MODEL_EXPLICIT implements the model interface and provides the function value and the first and second derivatives with respect to the parameters of an explicitly given model function.
%
% MODEL_EXPLICIT Methods:
%    GET_M - returns the result of the model function with parameter P and
%            experimental design T
%    GET_DP_M - returns the first derivative of the model function with
%               parameter P and experimental design T
%    GET_DPDP_M - returns the second derivative of the model function with
%                 parameter P and experimental design T
%
% see also MODEL
%

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
        f_sym;
        dp_f_sym;
        dpdp_f_sym;
        
        p_sym;
        t_sym;
    end
    
    methods (Access = public)
        
        function this = model_explicit(f, p, t)
        % MODEL_EXPLICIT creates a MODEL_EXPLICIT object.
        %
        % Example:
        %     OBJ = MODEL_EXPLICIT(F, P, T)
        %
        % Input:
        %     F: the explicit formula of the model function as string or 
        %        symbolic function. F may depend on P and T.
        %     P: the variables of the parameters P as string or symbolic vector
        %     T: the variables of the experimental design T as string or symbolic vector
        %
        % Output:
        %     OBJ: a MODEL_EXPLICIT object with the passed configurations
        %
        
            this.f_sym = util.make_sym(f);
            this.p_sym = util.make_sym(p);
            this.t_sym = util.make_sym(t);
            
            this.dp_f_sym = simplify(jacobian(this.f_sym, this.p_sym));
            this.dpdp_f_sym = simplify(jacobian(this.dp_f_sym, this.p_sym));
        end
        
        function M = get_M(this, p, t)        
        % GET_M returns the result of the model function with parameter P and experimental design T.
        %
        % Example:
        %     M = MODEL_EXPLICIT_OBJECT.GET_M(P, T)
        %
        % Input:
        %     P: the parameter values
        %     T: the experimental design values
        %
        % Output:
        %     M: the result of the model function with parameter P and experimental design T
        %
        
            M = this.substitute(this.f_sym, p, t);
        end
        
        function dp_M = get_dp_M(this, p, t)
        % GET_DP_M returns the first derivative of the model function with parameter P and experimental design T.
        %
        % Example:
        %     M = MODEL_EXPLICIT_OBJECT.GET_DP_M(P, T)
        %
        % Input:
        %     P: the parameter values
        %     T: the experimental design values
        %
        % Output:
        %     M: the first derivative of the model function with parameter P and experimental design T
        %
        
            dp_M = this.substitute(this.dp_f_sym, p, t);
        end
        
        function dpdp_M = get_dpdp_M(this, p, t)
        % GET_DP_M returns the second derivative of the model function with parameter P and experimental design T.
        %
        % Example:
        %     M = MODEL_EXPLICIT_OBJECT.GET_DPDP_M(P, T)
        %
        % Input:
        %     P: the parameter values
        %     T: the experimental design values
        %
        % Output:
        %     M: the second derivative of the model function with parameter P and experimental design T
        %
         
            dpdp_M = this.substitute(this.dpdp_f_sym, p, t);
        end
        
    end
    
    methods (Access = protected)
        
        function result = substitute(this, f_sym, p, t)
        % SUBSTITUTE substitutes the values P and T in F_SYM and returns the result.
        %
        % Example:
        %     RESULT = MODEL_EXPLICIT_OBJECT.SUBSTITUTE(F_SYM, P, T)
        %
        % Input:
        %     F_SYM: the symbolic formula
        %     P: the parameters value
        %     T: the experimental design values
        %
        % Output:
        %     RESULT: the result of the substitution of P and T in F_SYM
        %
                        
            p = util.make_row_vector(p);
            if not(all(size(p) == size(this.p_sym)))                
                error(util.get_message_identifier('model_explicit', 'substitute', 'wrong_size'), ['The vector p must have size ', mat2str(size(this.p_sym)), ' but its size is ', mat2str(size(p)), '.']);
            end
            if not(all(size(t) == size(this.t_sym)))                
                error(util.get_message_identifier('model_explicit', 'substitute', 'wrong_size'), ['The vector t must have size ', mat2str(size(this.t_sym)), ' but its size is ', mat2str(size(t)), '.']);
            end
            
            tmp_sym = subs(f_sym, this.p_sym, p);
            result = subs(tmp_sym, this.t_sym, t);
            result = double(result);
        end
        
    end
    
end
