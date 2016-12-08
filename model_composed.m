classdef model_composed < model
% MODEL_COMPOSED implements the model interface and provides the function value and the first and second derivatives with respect to the parameters of an explicitly given (symbolic) model function including an inner model.
%
% MODEL_COMPOSED Methods:
%    GET_M - returns the result of the model function with
%            model parameters P and experimental design X.
%    GET_DP_M - returns the first derivative of the model function with
%               model parameters P and experimental design X.
%    GET_DPDP_M - returns the second derivative of the model function with
%                 model parameters P and experimental design X.
%
% see also MODEL
%

%{
---------------------------------------------------------------------------
    Copyright (C) 2010-2016 Joscha Reimer jor@informatik.uni-kiel.de

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
    
    properties (Access = protected)
        p_sym;
        x_sym;
        
        f_sym;
        dp_f_sym;
        dpdp_f_sym;
        
        y_sym;
        dp_y_sym;
        dpdp_y_sym;
        
        inner_model;
    end
    
    methods (Access = public)
        
        function this = model_composed(f, p, x, inner_model_variable, inner_model_object)
        % MODEL_COMPOSED creates a MODEL_COMPOSED object.
        %
        % Example:
        %     OBJ = MODEL_COMPOSED(F, P, X)
        %
        % Input:
        %     F: the explicit formula of the model function as string or 
        %        symbolic function. F may depend on P and X.
        %     P: the variables of the model parameters as string or symbolic vector
        %     X: the variables of the experimental design as string or symbolic vector
        %     INNER_MODEL_VARIABLE: the variable used for the inner model value as string or symbolic vector
        %     INNER_MODEL_OBJECT: the inner model implementing the MODEL interface
        %
        % Output:
        %     OBJ: a MODEL_COMPOSED object with the passed configurations
        %
        
            % prepare p and t
            p = util.make_column_vector(util.make_sym(p));
            x = util.make_column_vector(util.make_sym(x));
            n = length(p);
            
            % prepare sym var for inner model and associated derivatives
            y = util.make_sym(inner_model_variable);
            dp_y = sym(['dp_' char(y)], [1, n]);
            dpdp_y = sym(['dpdp_' char(y)], [n, n]);
            
            % prepare f
            f = simplify(util.make_sym(f));
            
            % first derivative of f (regarding p)
            function d = der(f, x)
                d = simplify(jacobian(f, x));
            end
            
            dp_f_ = der(f, p);
            dy_f_ = der(f, y);
            
            dp_f = dp_f_ + dy_f_ * dp_y;
            
            % second derivative of f (regarding p)
            dpdp_f_ = der(dp_f_, p);
            dydp_f_ = der(dp_f_, y);
            dydy_f_ = der(dy_f_, y);
            dpdp_f = dpdp_f_;
            dpdp_f = dpdp_f + dydp_f_ * dp_y + dp_y' * dydp_f_';
            dpdp_f = dpdp_f + dp_y' * dydy_f_ * dp_y;
            dpdp_f = dpdp_f + dy_f_ * dpdp_y;
            
            % store needed values
            this.p_sym = p;
            this.x_sym = x;
            
            this.y_sym = y;
            this.dp_y_sym = dp_y;
            this.dpdp_y_sym = dpdp_y;
            
            this.f_sym = f;
            this.dp_f_sym = dp_f;
            this.dpdp_f_sym = dpdp_f;

            this.inner_model = inner_model_object;
        end
        
        function M = get_M(this, p, x)
        % GET_M returns the result of the model function with model parameters P and experimental design X.
        %
        % Example:
        %     M = MODEL_OBJECT.GET_M(P, X)
        %
        % Input:
        %     P: the model parameter values
        %     X: the experimental design values
        %
        % Output:
        %     M: the result of the model function with model parameters P and experimental design X
        %
        
            y = this.inner_model.get_M(p, x);
            
            M = subs(this.f_sym, this.y_sym, y);
            M = this.substitute(M, p, x);
            M = double(M);
        end
        
        function dp_M = get_dp_M(this, p, x)
        % GET_DP_M returns the first derivative of the model function with model parameters P and experimental design X.
        %
        % Example:
        %     M = MODEL_OBJECT.GET_DP_M(P, X)
        %
        % Input:
        %     P: the model parameter values
        %     X: the experimental design values
        %
        % Output:
        %     M: the first derivative of the model function with model parameters P and experimental design X
        %
        
            y = this.inner_model.get_M(p, x);
            dp_y = this.inner_model.get_dp_M(p, x);
            
            dp_M = this.dp_f_sym;
            dp_M = subs(dp_M, this.y_sym, y);
            dp_M = subs(dp_M, this.dp_y_sym, dp_y');
            dp_M = this.substitute(dp_M, p, x);
            dp_M = double(dp_M);
        end
        
        function dpdp_M = get_dpdp_M(this, p, x)
        % GET_DPDP_M returns the second derivative of the model function with model parameters P and experimental design X.
        %
        % Example:
        %     M = MODEL_OBJECT.GET_DPDP_M(P, X)
        %
        % Input:
        %     P: the model parameter values
        %     X: the experimental design values
        %
        % Output:
        %     M: the second derivative of the model function with model parameters P and experimental design X
        %
        
            y = this.inner_model.get_M(p, x);
            dp_y = this.inner_model.get_dp_M(p, x);
            dpdp_y = this.inner_model.get_dpdp_M(p, x);
            
            dpdp_M = this.dpdp_f_sym;
            dpdp_M = subs(dpdp_M, this.y_sym, y);
            dpdp_M = subs(dpdp_M, this.dp_y_sym, dp_y');
            dpdp_M = subs(dpdp_M, this.dpdp_y_sym, dpdp_y);
            dpdp_M = this.substitute(dpdp_M, p, x);
            dpdp_M = double(dpdp_M);
        end
        
    end
    
    methods (Access = protected)
        
        function result = substitute(this, f_sym, p, x)
        % SUBSTITUTE substitutes the values P and X in F_SYM and returns the result.
        %
        % Example:
        %     RESULT = MODEL_OBJECT.SUBSTITUTE(F_SYM, P, X)
        %
        % Input:
        %     F_SYM: the symbolic formula
        %     P: the model parameter values
        %     X: the experimental design values
        %
        % Output:
        %     RESULT: the result of the substitution of P and X in F_SYM
        %
                        
            p = util.make_column_vector(p);
            if not(all(size(p) == size(this.p_sym)))                
                error(util.get_message_identifier('model_explicit', 'substitute', 'wrong_size'), ['The vector p must have size ', mat2str(size(this.p_sym)), ' but its size is ', mat2str(size(p)), '.']);
            end
            x = util.make_column_vector(x);
            if not(all(size(x) == size(this.x_sym)))                
                error(util.get_message_identifier('model_explicit', 'substitute', 'wrong_size'), ['The vector x must have size ', mat2str(size(this.x_sym)), ' but its size is ', mat2str(size(x)), '.']);
            end
            
            tmp_sym = subs(f_sym, this.p_sym, p);
            result = subs(tmp_sym, this.x_sym, x);
        end
        
    end
    
end
