classdef model_fd < model
% MODEL_FD implements the model interface and provides the first and second derivatives with respect to the parameters by finite differences approximations.
%
% MODEL_FD Methods:
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

    properties (Access = protected)
        f = [];
        
        p = [];
        x = [];
        dp_M = [];
        dpdp_M = [];
    end
        
    methods (Access = public)
        
        function this = model_fd(f)
        % MODEL_FD creates a MODEL_FD object.
        %
        % Example:
        %     OBJ = MODEL_FD(F)
        %
        % Input:
        %     F: the function to vcalculate the model output
        %
        % Output:
        %     OBJ: a MODEL_FD object with the passed configurations
        %
        
            this.f = f;
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
        
            M = this.f(p, x);
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
        
             if ~ isequal(this.p, p) || ~ isequal(this.x, x) || isempty(this.dp_M)
                dp_M = util.approximate_Jacobian(@(p) (this.get_M(p, x)), p);
                                
                this.dp_M = dp_M;
                this.p = p;
                this.x = x;
            else
                dp_M = this.dp_M;
            end           
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
         
            if ~ isequal(this.p, p) || ~ isequal(this.x, x) || isempty(this.dpdp_M)
                dp_M = this.get_dp_M(p, x);
                dpdp_M = util.approximate_Hessian(@(p) (this.get_M(p, x)), p, dp_M);
                
                this.dp_M = dp_M;
                this.dpdp_M = dpdp_M;
                this.p = p;
                this.x = x;
            else
                dpdp_M = this.dpdp_M;
            end
        end
        
    end
    
end

