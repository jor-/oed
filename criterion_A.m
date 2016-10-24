classdef criterion_A < criterion
% CRITERION_A is the quality criterion that uses the average variance as the quality.
% In other words it holds PHI = 1/n * trace(C) for the quality PHI and
% the covariance matrix C if n x n is the dimension of C.
%
% CRITERION_A Methods:
%    GET_PHI - returns the quality of the covariance matrix C.
%    GET_DCD_PHI - returns the derivation in the direction D of the
%                  quality of C.
%    GET_DCD_DCD_PHI - returns the derivation with respect to C in the
%                      direction E of the derivation in the direction D
%                      of the quality of C.
%    GET_DDD_DCD_PHI - returns the derivation with respect to D in the
%                      direction E of the derivation in the direction D
%                      of the quality of C.
%
% see also CRITERION
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
    
    methods (Static)
        function phi = get_phi(C)
        % GET_PHI returns the quality of the covariance matrix C.
        %
        % Example:
        %     PHI = CRITERION_OBJECT.GET_PHI(C)
        %
        % Input:
        %     C: the covariance matrix C
        %
        % Output:
        %     PHI: the quality of the covariance matrix C
        %
            size_C = size(C);
            phi = trace(C) / size_C(1);
        end
        
        function dCd_phi = get_dCd_phi(C, D)
        % GET_DPD_PHI returns the derivation in the direction D of the quality of C.
        %
        % Example:
        %     DCD_PHI = CRITERION_OBJECT.GET_DCD_PHI(C, D)
        %
        % Input:
        %     C: the covariance matrix C
        %     D: the direction D of the derivate
        %
        % Output:
        %     DCD_PHI: the derivation in the direction D of the quality of C
        %
            size_C = size(C);
            dCd_phi = trace(D) / size_C(1);
        end
        
        function dCd_dCd_phi = get_dCd_dCd_phi(C, D, E)
        % GET_DCD_DCD_PHI returns the derivation with respect to C in the direction E of the derivation in the direction D of the quality of C.
        %
        % Example:
        %     DCD_DCD_PHI = CRITERION_OBJECT.GET_DWD_PHI(C, D, E)
        %
        % Input:
        %     C: the covariance matrix C
        %     D: the direction D of the first derivate
        %     E: the direction E of the second derivate
        %
        % Output:
        %     DCD_DCD_PHI: the derivation with respect to C in the direction
        %                  E of the derivation in the direction D of the
        %                  quality of C 
        %
            dCd_dCd_phi = 0;
        end
        
        function dDd_dCd_phi = get_dDd_dCd_phi(C, D, E)
        % GET_DDD_DCD_PHI returns the derivation with respect to D in the direction E of the derivation in the direction D of the quality of C.
        %
        % Example:
        %     DDD_DCD_PHI = CRITERION_OBJECT.GET_DWD_PHI(C, D, E)
        %
        % Input:
        %     C: the covariance matrix C
        %     D: the direction D of the first derivate
        %     E: the direction E of the second derivate
        %
        % Output:
        %     DDD_DCD_PHI: the derivation with respect to D in the direction
        %                  E of the derivation in the direction D of the
        %                  quality of C 
        %
            size_C = size(C);
            dDd_dCd_phi = trace(E) / size_C(1);            
        end
    end
    
end

