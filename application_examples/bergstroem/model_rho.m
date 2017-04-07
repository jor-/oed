classdef model_rho < model_ivp
% MODEL_RHO implements the model interface and provides the function value and the first and second derivatives with respect to the model parameters.
%
% MODEL_RHO Methods:
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
        
    methods (Access = public)
        
        function this = model_rho(R)
        % MODEL_RHO creates a MODEL_RHO object.
        %
        % Example:
        %     OBJ = MODEL_RHO(R)
        %
        % Input:
        %     R: the model constant R
        %
        % Output:
        %     OBJ: a MODEL_RHO object with the passed configurations
        %

            % model parameters
            rho0    = sym('rho0'); 
            m0      = sym('m0');
            U0      = sym('U0'); 
            Omega0  = sym('Omega0'); 
            C       = sym('C'); 
            Qv      = sym('Qv');
            p = [rho0; m0; U0; Omega0; C; Qv];

            % differentiation variable
            epsilon       = sym('epsilon');
            epsilon_interval = [0; 1];

            % independent model variables (without differentiation variable)
            T             = sym('T');
            epsilon_dot   = sym('epsilon_dot'); 
            x = [T; epsilon_dot;];

            % right side of differential equation
            rho = sym('rho');
            Omega = Omega0 + C * exp(-m0*Qv/R/(T+273)) * epsilon_dot^(-m0);
            f = U0 * rho^(1/2) - Omega * rho;

            % inital value problem model for rho
            this = this@model_ivp(f, p, rho, rho0, epsilon, epsilon_interval, x);
        end
        
    end
    
end
