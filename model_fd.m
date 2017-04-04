classdef model_fd < model
% MODEL_FD implements the model interface and provides the first and second derivatives with respect to the parameters by finite differences approximations.
%
% MODEL_FD Methods:
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
        p = [];
        t = [];
        dp_M = [];
        dpdp_M = [];
    end
        
    methods (Access = public)
        
        function dp_M = get_dp_M(this, p, t)
        % GET_DP_M returns the first derivative of the model function with parameter P and experimental design T.
        %
        % Example:
        %     M = MODEL_FD_OBJECT.GET_DP_M(P, T)
        %
        % Input:
        %     P: the parameter values
        %     T: the experimental design values
        %
        % Output:
        %     M: the first derivative of the model function with parameter P and experimental design T
        %
        
             if ~ isequal(this.p, p) || ~ isequal(this.t, t) || isempty(this.dp_M)
                dp_M = util.approximate_Jacobian(@(p) (this.get_M(p, t)), p);
                                
                this.dp_M = dp_M;
                this.p = p;
                this.t = t;
            else
                dp_M = this.dp_M;
            end           
        end
        
        function dpdp_M = get_dpdp_M(this, p, t)
        % GET_DP_M returns the second derivative of the model function with parameter P and experimental design T.
        %
        % Example:
        %     M = MODEL_FD_OBJECT.GET_DPDP_M(P, T)
        %
        % Input:
        %     P: the parameter values
        %     T: the experimental design values
        %
        % Output:
        %     M: the second derivative of the model function with parameter P and experimental design T
        %
         
            if ~ isequal(this.p, p) || ~ isequal(this.t, t) || isempty(this.dpdp_M)
                dp_M = this.get_dp_M(p, t);
                dpdp_M = util.approximate_Hessian(@(p) (this.get_M(p, t)), p, dp_M);
                
                this.dp_M = dp_M;
                this.dpdp_M = dpdp_M;
                this.p = p;
                this.t = t;
            else
                dpdp_M = this.dpdp_M;
            end
        end
        
    end
    
end

