classdef model_C3 < model_C2
% MODEL_C3 models the concentration C of the suspended sediment above the marsh surface with three parameters.
%
% MODEL_C3 Methods:
%   GET_T_SPAN - returns the span in which the marsh is flooded.
%   GET_M - returns the concentration C of the suspended sediment.
%   GET_DP_M - returns the derivation of the concentration C of the suspended sediment with respect to the parameters.
%   GET_DPDP_M - returns the second derivation of the concentration C of the suspended sediment with respect to the parameters.
%   SET_E - sets the elevation of the marsh surface.
%   SET_A - sets the factor A which influences the flooding of the marsh surface.
%   SET_B - sets the factor B which influences the flooding of the marsh surface.
%   SET_X0 - sets the factor X0 which influences the flooding of the marsh surface.
%   SET_hMHW - sets the mean high-water level.
%   SET_hHW - sets the the certain high-water level.
%   SET_DEBUG - enables or disables the debug output. 
%
% see also MODEL_C2
%

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
    
    methods (Access = public)
        
        % ******************* CONSTRUCTOR ******************* %
        function this = model_C3(E, a, b, x0, hMHW, hHW)
        % MODEL_C3 creates a MODEL_C3 object.
        %
        % Example:
        %     OBJ = MODEL_C3(E, a, b, x0, hMHW, hHW)
        %
        % Input:
        %     E: the elevation of the marsh surface
        %        (optional, can later be set by SET_E)
        %     a: a factor which influences the flooding of the marsh surface
        %        (optional, can later be set by SET_A)
        %     b: a factor which influences the flooding of the marsh surface
        %        (optional, can later be set by SET_B)
        %     x0: a factor which influences the flooding of the marsh surface
        %         (optional, can later be set by SET_X0)
        %     hMHW: the the mean high-water level
        %           (optional, can later be set by SET_hMHW)
        %     hHW: the certain high-water level
        %          (optional, can later be set by SET_hHW)
        %
        % Output:
        %     OBJ: a MODEL_C3 object with the passed configurations
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hHMW and SET_hHW
        %
        
            if nargin == 0
                E = [];
            end
            if nargin <= 1
                a = [];
            end
            if nargin <= 2
                b = [];
            end
            if nargin <= 3
                x0 = [];
            end
            if nargin <= 4
                hMHW = [];
            end
            if nargin <= 5
                hHW = [];
            end
        
            this = this@model_C2(E, a, b, x0, hMHW, hHW);
        end        
        
    end
    
    
    methods (Access = protected)
        
        function p_sym = get_p_sym(this)
        % GET_P_SYM returns the symbolic formula for the parameters P.
        %
        % Example:
        %     P_SYM = MODEL_C3.GET_P_SYM()
        %
        % Output:
        %     P_SYM: the symbolic formula for the parameters P
        %
        
            p_str = {'k'; 'q'; 'r'};

            for i=1:length(p_str)
                p_sym(i) = sym(p_str{i});
            end
        end    
        
        function f_sym = get_f_sym(this)
        % GET_F_SYM returns the symbolic formula for the derivation of the concentration C with respect to the point in time T.
        %
        % Example:
        %     F_SYM = MODEL_C3.GET_F_SYM()
        %
        % Output:
        %     F_SYM: the symbolic formula for the derivation of the concentration C with respect to the point in time T
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % the corresponding SET methods or the input.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hHMW and SET_hHW
        %
        
            t_sym = this.get_t_sym();
            h_sym = this.get_h_sym();            

            C_sym = this.get_C_sym;            
            p_sym = this.get_p_sym();

            hHW = this.get_hHW();
            E = this.get_E();

            C0_sym = p_sym(1) * (hHW - E);
            ws_sym = p_sym(2) * C0_sym^p_sym(3);
            
            dt_h_sym = simplify(jacobian(h_sym, t_sym));
            dt_hp_sym = simplify((dt_h_sym + abs(dt_h_sym)) / 2);

            f_sym = simplify((-ws_sym * C_sym + dt_hp_sym * (C0_sym-C_sym)) / (h_sym - E));
        end  
                
        function C0_sym = get_C0_sym(this)
        % GET_C0_sym returns the symbolic formula for the initial value of the inital value problem.
        %
        % Example:
        %     C0_sym = MODEL_C3.GET_C0_sym()
        %
        % Output:
        %     C0_sym: the symbolic formula for the initial value of the
        %             inital value problem
        %
        
            p_sym = this.get_p_sym();

            hHW = this.get_hHW();
            E = this.get_E();

            C0_sym = p_sym(1) * (hHW - E);
        end
        
    end
    
end


%#ok<*PROP>
