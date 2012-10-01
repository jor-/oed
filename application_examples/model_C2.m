classdef model_C2 < model
% MODEL_C2 models the concentration C of the suspended sediment above the marsh surface with two parameters.
%
% MODEL_C2 Methods:
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
% see also MODEL
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
    
    properties (Access = private)
        E;
        a;
        b;
        x0;
        hMHW;
        hHW;
        
        p;
        t;
        
        t_span;        
        ivp;
        
        debug;
    end
    
    events (ListenAccess = protected, NotifyAccess = protected)
        event_E_changed;    %is triggered when E was changed
        event_a_changed;    %is triggered when a was changed
        event_b_changed;    %is triggered when b was changed
        event_x0_changed;   %is triggered when x0 was changed
        event_hMHW_changed; %is triggered when hMHW was changed
        event_hHW_changed;  %is triggered when hHW was changed
        event_p_changed;    %is triggered when p was changed
        event_t_changed;    %is triggered when t was changed
    end
    
    methods (Access = public)
        
        % ******************* CONSTRUCTOR ******************* %
        
        function this = model_C2(E, a, b, x0, hMHW, hHW)
        % MODEL_C2 creates a MODEL_C2 object.
        %
        % Example:
        %     OBJ = MODEL_C2(E, a, b, x0, hMHW, hHW)
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
        %     OBJ: a MODEL_C2 object with the passed configurations
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            if nargin >= 1
                this.set_E(E);
            end
            if nargin >= 2
                this.set_a(a);
            end
            if nargin >= 3
                this.set_b(b);
            end
            if nargin >= 4
                this.set_x0(x0);
            end
            if nargin >= 5
                this.set_hMHW(hMHW);
            end
            if nargin >= 6
                this.set_hHW(hHW);
            end            
                    
            this.set_debug(0);
            
            
            addlistener(this, 'event_E_changed', @(src, evnt_data) remove_calculations('E'));
            addlistener(this, 'event_a_changed', @(src, evnt_data) remove_calculations('a'));
            addlistener(this, 'event_b_changed', @(src, evnt_data) remove_calculations('b'));
            addlistener(this, 'event_x0_changed', @(src, evnt_data) remove_calculations('x0'));            
            addlistener(this, 'event_hMHW_changed', @(src, evnt_data) remove_calculations('hMHW'));
            addlistener(this, 'event_hHW_changed', @(src, evnt_data) remove_calculations('hHW'));
            addlistener(this, 'event_p_changed', @(src, evnt_data) remove_calculations('p'));
            addlistener(this, 'event_t_changed', @(src, evnt_data) remove_calculations('t'));
            
            function remove_calculations(changed)
                if strfind('E a b x0 hHMW hHW',  changed)
                    this.t_span = [];
                    this.ivp = [];
                end
                if strfind('p',  changed)
                    this.ivp = [];
                end
            end
        end        
        
        
        % ******************* T SPAN ******************* %
        
        function t_span = get_t_span(this, hHW)
        % GET_T_SPAN returns the span in which the marsh is flooded.
        %
        % Example:
        %     T_SPAN = MODEL_C2.GET_T_SPAN(hHW)
        %
        % Input:
        %     hHW: the certain high-water level
        %          (optional, if hHW has been set previously)
        %
        % Output:
        %     T_SPAN: a vector specifying the span in which the marsh is
        %             flooded. The first component is the initial point and
        %             the second component the end point.
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % the corresponding SET methods or the input.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            if nargin >= 2
                this.set_hHW(hHW);
            end
            
            if isempty(this.t_span)
                h_sym = this.get_h_sym();
                E = this.get_E();

                % calculate flooded time span
                r = subs(solve(h_sym - E));

                if length(r) ~= 2 || not(isreal(r))
                    error('Es konnten nicht die Zeitspanne bestimmt werden, in der die Wiese ueberflutet ist. Vielleicht ist die Wiese immer oder nie ueberflutet? Stellen Sie sicher, dass "a > E - h_HW + h_MHW" gilt.');
                end

                r = sort(r);
                t_span = [r(1) r(2)] + [1 -1];
                this.t_span = t_span;
            else
                t_span = this.t_span;
            end
        end        
                
        
        % ******************* MODEL Output ******************* %
        
        function M = get_M(this, p, t)
        % GET_M returns the concentration C of the suspended sediment.
        %
        % Example:
        %     M = MODEL_C2.GET_M(P, T)
        %
        % Input:
        %     P: the parameters
        %     T: the measuring point
        %
        % Output:
        %     M: the concentration C of the suspended sediment for the passed parameters P at the measuring point T
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % or the corresponding SET methods.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            if nargin >= 2
                this.set_p(p);
            end
            if nargin >= 3
                this.set_t(t);
            end
            
            ivp = this.get_ivp();
            p = this.get_p();
            t = this.get_t();            
            
            M = ivp.get_M(p, t);
        end        
        
        function dp_M = get_dp_M(this, p, t)
        % GET_DP_M returns the derivation of the concentration C of the suspended sediment with respect to the parameters.
        %
        % Example:
        %     M = MODEL_C2.GET_DP_M(P, T)
        %
        % Input:
        %     P: the parameters
        %     T: the measuring point
        %
        % Output:
        %     M: the derivation of the concentration C of the suspended sediment with respect to the parameters
        %        for the passed parameters P at the measuring point T
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % or the corresponding SET methods.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            if nargin >= 2
                this.set_p(p);
            end
            if nargin >= 3
                this.set_t(t);
            end
            
            ivp = this.get_ivp();
            p = this.get_p();
            t = this.get_t();
            
            dp_M = ivp.get_dp_M(p, t);
        end
                
        function dpdp_M = get_dpdp_M(this, p, t)
        % GET_DPDP_M returns the second derivation of the concentration C of the suspended sediment with respect to the parameters.
        %
        % Example:
        %     M = MODEL_C2.GET_DPDP_M(P, T)
        %
        % Input:
        %     P: the parameters
        %     T: the measuring point
        %
        % Output:
        %     M: the second derivation of the concentration C of the suspended sediment with respect to the parameters
        %        for the passed parameters P at the measuring point T
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % or the corresponding SET methods.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
        
            if nargin >= 2
                this.set_p(p);
            end
            if nargin >= 3
                this.set_t(t);
            end
            
            ivp = this.get_ivp();
            p = this.get_p();
            t = this.get_t();
            
            dpdp_M = ivp.get_dpdp_M(p, t);
        end
                
                
        % ******************* THE SETTER FOR THE INPUTS ******************* %
        
        function set_E(this, E)
        % SET_E sets the elevation of the marsh surface.
        %
        % Example:
        %     MODEL_C2.SET_E(E)
        %
        % Input:
        %     E: the elevation of the marsh surface
        %
        % see also GET_T_SPAN, GET_M, GET_DP_M and GET_DPDP_M
        %
        
            if not(isequal(this.E, E))
                this.E = E;
                notify(this, 'event_E_changed');
            end
        end
        
        function set_a(this, a)
        % SET_A sets the factor A which influences the flooding of the marsh surface.
        %
        % Example:
        %     MODEL_C2.SET_A(A)
        %
        % Input:
        %     A: the factor A which influences the flooding of the marsh surface
        %
        % see also GET_T_SPAN, GET_M, GET_DP_M and GET_DPDP_M
        %
        
            if not(isequal(this.a, a))
                this.a = a;
                notify(this, 'event_a_changed');
            end
        end
        
        function set_b(this, b)
        % SET_B sets the factor B which influences the flooding of the marsh surface.
        %
        % Example:
        %     MODEL_C2.SET_B(B)
        %
        % Input:
        %     B: the factor B which influences the flooding of the marsh surface
        %
        % see also GET_T_SPAN, GET_M, GET_DP_M and GET_DPDP_M
        %
        
            if not(isequal(this.b, b))
                this.b = b;
                notify(this, 'event_b_changed');
            end
        end
        
        function set_x0(this, x0)
        % SET_X0 sets the factor X0 which influences the flooding of the marsh surface.
        %
        % Example:
        %     MODEL_C2.SET_X0(X0)
        %
        % Input:
        %     X0: the factor X0 which influences the flooding of the marsh surface
        %
        % see also GET_T_SPAN, GET_M, GET_DP_M and GET_DPDP_M
        %
        
            if not(isequal(this.x0, x0))
                this.x0 = x0;
                notify(this, 'event_x0_changed');
            end
        end
        
        function set_hMHW(this, hMHW)
        % SET_hMHW sets the mean high-water level.
        %
        % Example:
        %     MODEL_C2.SET_hMHW(hMHW)
        %
        % Input:
        %     hMHW: the mean high-water level
        %
        % see also GET_T_SPAN, GET_M, GET_DP_M and GET_DPDP_M
        %
        
            if not(isequal(this.hMHW, hMHW))
                this.hMHW = hMHW;
                notify(this, 'event_hMHW_changed');
            end
        end
        
        function set_hHW(this, hHW)
        % SET_hHW sets the the certain high-water level.
        %
        % Example:
        %     MODEL_C2.SET_hHW(hHW)
        %
        % Input:
        %     hHW: the the certain high-water level
        %
        % see also GET_T_SPAN, GET_M, GET_DP_M and GET_DPDP_M
        %
        
            if not(isequal(this.hHW, hHW))
                this.hHW = hHW;
                notify(this, 'event_hHW_changed');
            end
        end    
        
        function set_debug(this, debug)
        % SET_DEBUG enables or disables the debug output.
        %
        % Example:
        %     MODEL_C2.SET_DEBUG(DEBUG)
        %
        % Input:
        %     DEBUG: zero for disabling the debug output
        %            a positive integer for enabling the debug output
        %
        
            if not(isempty(this.ivp))
                this.ivp.set_debug(debug);
            end
            this.debug = debug;
        end
        
    end
    
    
    
    methods (Access = protected)
        
        % ******************* INTERNAL GETTER AND SETTER ******************* %
        
        function t_sym = get_t_sym(this)
        % GET_T_SYM returns the symbolic formula for the point in time T.
        %
        % Example:
        %     T_SYM = MODEL_C2.GET_T_SYM()
        %
        % Output:
        %     T_SYM: the symbolic formula for the point in time
        %
        
            t_sym = sym('t');
        end
        
        function C_sym = get_C_sym(this)
        % GET_C_SYM returns the symbolic formula for the concentration C.
        %
        % Example:
        %     T_SYM = MODEL_C2.GET_T_SYM()
        %
        % Output:
        %     T_SYM: the symbolic formula for the concentration C
        %
        
            C_sym = sym('C');    
        end
                
        function p_sym = get_p_sym(this)
        % GET_P_SYM returns the symbolic formula for the parameters P.
        %
        % Example:
        %     P_SYM = MODEL_C2.GET_P_SYM()
        %
        % Output:
        %     P_SYM: the symbolic formula for the parameters P
        %
        
            p_str = {'C0'; 'ws'};

            for i=1:length(p_str)
                p_sym(i) = sym(p_str{i});
            end
        end 
        
        function h_sym = get_h_sym(this)
        % GET_H_SYM returns the symbolic formula for the water surface elevation.
        %
        % Example:
        %     H_SYM = MODEL_C2.GET_H_SYM()
        %
        % Output:
        %     H_SYM: the symbolic formula for the water surface elevation
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % the corresponding SET methods or the input.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            t_sym = this.get_t_sym();
            a = this.get_a();
            b= this.get_b();
            x0 = this.get_x0();            
            hMHW = this.get_hMHW();
            hHW = this.get_hHW();
            h_sym = a / (1 + ((t_sym - x0) / b).^2) + hHW - hMHW;
        end
       
        function f_sym = get_f_sym(this)
        % GET_F_SYM returns the symbolic formula for the derivation of the concentration C with respect to the point in time T.
        %
        % Example:
        %     F_SYM = MODEL_C2.GET_F_SYM()
        %
        % Output:
        %     F_SYM: the symbolic formula for the derivation of the concentration C with respect to the point in time T
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % the corresponding SET methods or the input.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            t_sym = this.get_t_sym();
            h_sym = this.get_h_sym();            

            C_sym = this.get_C_sym;            
            p_sym = this.get_p_sym();

            E = this.get_E();

            C0_sym = p_sym(1);
            ws_sym = p_sym(2);
            
            dt_h_sym = simplify(jacobian(h_sym, t_sym));
            dt_hp_sym = simplify((dt_h_sym + abs(dt_h_sym)) / 2);

            f_sym = simplify((-ws_sym * C_sym + dt_hp_sym * (C0_sym-C_sym)) / (h_sym - E));
        end
        
        function C0_sym = get_C0_sym(this)
        % GET_C0_sym returns the symbolic formula for the initial value of the initial value problem.
        %
        % Example:
        %     C0_sym = MODEL_C2.GET_C0_sym()
        %
        % Output:
        %     C0_sym: the symbolic formula for the initial value of the
        %             initial value problem
        %
        
            p_sym = this.get_p_sym();
            
            C0_sym = p_sym(1);
        end
        
        
        
        function ivp = get_ivp(this)
        % GET_IVPDS returns the IVP_DERIVATION_SOLVER-object which provides the concentration C and its first and second derivation with respect to the parameters P.
        %
        % Example:
        %     IVPDS = MODEL_C2.GET_IVPDS()
        %
        % Output:
        %     IVPDS: the IVP_DERIVATION_SOLVER-object which provides the concentration C 
        %            and its first and second derivation with respect to the parameters P
        %
        % E, A, B, X0, hHMW and hHW must be set before via the constructor,
        % the corresponding SET methods or the input.
        %
        % see also SET_E, SET_A, SET_B, SET_X0, SET_hMHW and SET_hHW
        %
        
            if isempty(this.ivp)
                t_sym = this.get_t_sym();
                C_sym = this.get_C_sym();            
                p_sym = this.get_p_sym();
                f_sym = this.get_f_sym();                
                t_span = this.get_t_span();
                y0_sym = this.get_C0_sym;
                
                ivp = model_ivp(t_sym, C_sym, p_sym, f_sym, t_span, y0_sym);
                ivp.set_debug(this.get_debug());
                
                this.ivp = ivp;
            else
                ivp = this.ivp;
            end
        end
        
        
               
        function E = get_E(this)
        % GET_E returns the elevation of the marsh surface.
        %
        % Example:
        %     E = MODEL_C2.GET_E()
        %
        % Output:
        %     E: the elevation of the marsh surface
        %
        % see also SET_E
        %
        
            if isempty(this.E)
                error(this.get_message_identifier('get_E', 'not_set'), 'E is not set.');
            else
                E = this.E;
            end 
        end        
        
        function a = get_a(this)
        % GET_A returns the factor A which influences the flooding of the marsh surface.
        %
        % Example:
        %     A = MODEL_C2.GET_A()
        %
        % Output:
        %     A: the factor A which influences the flooding of the marsh surface
        %
        % see also SET_A
        %
        
            if isempty(this.a)
                error(this.get_message_identifier('get_a', 'not_set'), 'a is not set.');
            else
                a = this.a;
            end 
        end
        
        function b = get_b(this)
        % GET_B returns the factor B which influences the flooding of the marsh surface.
        %
        % Example:
        %     B = MODEL_C2.GET_B()
        %
        % Output:
        %     B: the factor B which influences the flooding of the marsh surface
        %
        % see also SET_B
        %
        
            if isempty(this.b)
                error(this.get_message_identifier('get_b', 'not_set'), 'b is not set.');
            else
                b = this.b;
            end 
        end
        
        function x0 = get_x0(this)
        % GET_X0 returns the factor X0 which influences the flooding of the marsh surface.
        %
        % Example:
        %     X0 = MODEL_C2.GET_X0()
        %
        % Output:
        %     X0: the factor X0 which influences the flooding of the marsh surface
        %
        % see also SET_X0
        %
        
            if isempty(this.x0)
                error(this.get_message_identifier('get_x0', 'not_set'), 'x0 is not set.');
            else
                x0 = this.x0;
            end 
        end
        
        function hMHW = get_hMHW(this)
        % GET_hMHW returns the mean high-water level.
        %
        % Example:
        %     hMHW = MODEL_C2.GET_hMHW()
        %
        % Output:
        %     hMHW: the mean high-water level
        %
        % see also SET_hMHW
        %
        
            if isempty(this.hMHW)
                error(this.get_message_identifier('get_hMHW', 'not_set'), 'hMHW is not set.');
            else
                hMHW = this.hMHW;
            end 
        end
        
        function hHW = get_hHW(this)
        % GET_hHW returns the certain high-water level.
        %
        % Example:
        %     hHW = MODEL_C2.GET_hHW()
        %
        % Output:
        %     hHW: the certain high-water level
        %
        % see also SET_hHW
        %
        
            if isempty(this.hHW)
                error(this.get_message_identifier('get_hHW', 'not_set'), 'hHW is not set.');
            else
                hHW = this.hHW;
            end 
        end
        
        
        function p = get_p(this)
        % GET_P returns the parameters P.
        %
        % Example:
        %     P = MODEL_C2.GET_P()
        %
        % Output:
        %     P: the parameters P
        %
        % see also SET_P
        %
        
            if isempty(this.p)
                error(this.get_message_identifier('get_p', 'not_set'), 'p is not set.');
            else
                p = this.p;
            end 
        end
        
        function set_p(this, p)
        % SET_P sets the parameters P.
        %
        % Example:
        %     MODEL_C2.SET_P(P)
        %
        % Input:
        %     P: the parameters P
        %
        % see also GET_P
        %
        
            p = util.make_column_vector(p);
            
            if not(isequal(this.p, p))
                this.p = p;
                notify(this, 'event_p_changed');
            end
        end
        
        
        function t = get_t(this)
        % GET_P returns the point in time T.
        %
        % Example:
        %     T = MODEL_C2.GET_T()
        %
        % Output:
        %     T: the point in time T
        %
        % see also SET_T
        %
        
            if isempty(this.t)
                error(this.get_message_identifier('get_t', 'not_set'), 't is not set.');
            else
                t = this.t;
            end 
        end
        
        function set_t(this, t)
        % SET_T sets the measuring point T.
        %
        % Example:
        %     MODEL_C2.SET_T(T)
        %
        % Input:
        %     T: the measuring point T. The first component is the point in
        %        time. The optional second component is the certain
        %        high-water level hHW.
        %
        % see also GET_T and GET_hHW
        %
            
            if iscell(t)
               t = t{1}; 
            end
        
            if not(isequal(this.t, t(1)))
                this.t = t(1);
                notify(this, 'event_t_changed');
            end
            
            if length(t) == 2
                this.set_hHW(t(2));
            end
        end       
        
        
        function debug = get_debug(this)
        % GET_DEBUG returns the current debug-level.
        %
        % Example:
        %     DEBUG = MODEL_C2.GET_DEBUG()
        %
        % Output:
        %     DEBUG: the current debug-level
        %
        % see also SET_DEBUG
        %        
        
            debug = this.debug;
        end
        
    end
    
    methods (Access = protected, Static)
        
        function id = get_message_identifier(method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in methods of these object.
        %
        % Example:
        %     ID = MODEL_C2.GET_MESSAGE_IDENTIFIER(METHOD, MNEMONIC)
        %
        % Input:
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
        
            id = util.get_message_identifier('model_C2', method, mnemonic);
        end
        
    end
    
end


%#ok<*PROP>
