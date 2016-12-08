classdef model_bergstroem < model_composed
% MODEL_BERGSTROEM implements the model interface and provides the function value and the first and second derivatives with respect to the model parameters.
%
% MODEL_BERGSTROEM Methods:
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
        
        function this = model_bergstroem(c)
        % MODEL_BERGSTROEM creates a MODEL_BERGSTROEM object.
        %
        % Example:
        %     OBJ = MODEL_BERGSTROEM(F, P, T)
        %
        % Input:
        %     C: the model constants as a vector
        %        [alpha, b, R, c_T, T1, Gref]
        %
        % Output:
        %     OBJ: a MODEL_BERGSTROEM object with the passed configurations
        %

            % model constants
            alpha = c(1);
            b     = c(2);
            R     = c(3);
            C_T   = c(4);
            T1    = c(5);
            Gref  = c(6);

            % model parameters
            rho0    = sym('rho0'); 
            m0      = sym('m0');
            U0      = sym('U0'); 
            Omega0  = sym('Omega0'); 
            C       = sym('C'); 
            Qv      = sym('Qv');
            p = [rho0; m0; U0; Omega0; C; Qv];

            % independent model variables
            T             = sym('T'); 
            epsilon       = sym('epsilon');
            epsilon_dot   = sym('epsilon_dot'); 
            x = [T; epsilon_dot; epsilon];

            % inner model
            rho = sym('rho');
            rho_model = model_rho(R);
                       
            % function for sigma_f
            sigma_0 = alpha * Gref * b * sqrt(rho0);
            sigma_f = (1 - C_T * exp(-T1/T)) * (sigma_0 + alpha * Gref * b * sqrt(rho));
            
            % use model_composed constructor
            this = this@model_composed(sigma_f, p, x, rho, rho_model);
        end
        
    end
    
end
