classdef model < handle
% MODEL represents an interface for a model.
%
% MODEL Methods:
%    GET_M - returns the result of the model function
%    GET_DP_M - returns the first derivative of the model function 
%               with respect to the parameter p at point t	
%    GET_DPDP_M - returns the second derivative of the model function 
%                 with respect to the parameter p at point t

%{
---------------------------------------------------------------------------
    Copyright (C) 2010-2012 Joscha Reimer jor@informatik.uni-kiel.de

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
	
	methods (Abstract)
		
        M = get_M(this, p, t)
        % GET_M returns the result of the model function.
        %
        % Example:
        %     M = MODEL_OBJECT.GET_M(P, T)
        %
        % Input:
        %     P: the parameters
        %     T: the measurement
        %
        % Output:
        %     M: the result of the model function
        %
        
		dp_M = get_dp_M(this, p, t)
        % GET_DP_M returns the first derivative of the model function with respect to the parameters P at point T.
        %
        % Example:
        %     M = MODEL_OBJECT.GET_DP_M(P, T)
        %
        % Input:
        %     P: the parameters
        %     T: the measurement
        %
        % Output:
        %     M: the first derivative of the model function with respect to the parameters P at point T
        %
		
		dpdp_M = get_dpdp_M(this, p, t)
		% GET_DP_M returns the second derivative of the model function with respect to the parameters P at point T.
        %
        % Example:
        %     M = MODEL_OBJECT.GET_DPDP_M(P, T)
        %
        % Input:
        %     P: the parameters
        %     T: the measurement
        %
        % Output:
        %     M: the second derivative of the model function with respect to the parameters P at point T
        %
		
	end
	
end

