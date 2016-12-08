classdef model_ivp < model
% MODEL_IVP implements the model interface and provides the solution of an initial value problem and its first and second derivatives with respect to the parameters.
%
% MODEL_IVP Methods:
%   GET_M - returns the solution of the initial value problem with
%           parameter P.
%   GET_DP_M - returns the first derivative with respect to the
%              parameters of the solution of the initial value problem
%              with parameter P. 
%   GET_DPDP_M - returns the second derivative with respect to the
%                parameters of the solution of the initial value problem
%                with parameter P.
%   SET_DEBUG - enables or disables the debug output.
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
        p_sym;
        x_sym;
        t_sym;
        
        p;
        x;
        t;
        t_interval;
        
        f_sym;
        y_sym;
        y0_sym;
        y_p_x;
        y_p_x_t;
        
        dp_f_sym;
        dp_y_sym;
        dp_y0_sym;    
        dp_y_p_x;
        dp_y_p_x_t;
        
        dpdp_f_sym;
        dpdp_y_sym;
        dpdp_y0_sym;
        dpdp_y_p_x;
        dpdp_y_p_x_t;
        
        ode_solver_name;
        debug;
    end
    
    events (ListenAccess = protected, NotifyAccess = protected)
        event_p_changed;    %is triggered when p was changed
        event_t_changed;    %is triggered when t was changed
        event_x_changed;    %is triggered when x was changed
    end
    
    methods (Access = public)
        
        function this = model_ivp(f, p, y, y0, t, t_interval, x, ode_solver_name)
        % MODEL_IVP creates an MODEL_IVP object.
        %
        % Example:
        %     OBJ = MODEL_IVP(F, P, Y, Y0, T, T_INTERVAL, X)
        %
        % Input:
        %     F: the formula f of the derivative of y with respect to t as
        %        string or symbolic formula. f can depend on t, x, y and p.
        %        It holds (dy/dt)(t,x,p) = f(y,t,x,p).
        %     P: the variables of the parameters p as a cell array of
        %        strings or as symbolic vector
        %     Y: the dependent variable y as string or symbolic variable
        %     Y_0: the value of y at t_0 as string or symbolic variable.
        %          Y_0_SYM can depend on x and the parameters p. It holds
        %          y(t_0,x,p) = y_0(x,p).
        %     T: the independent variable t of the ordinary differential 
        %        equation as string or symbolic variable
        %     T_INTERVAL: a vector specifying the interval of integration,
        %             T_INTERVAL = [t_0, t_f]. The solver imposes the initial
        %             conditions at t0 and integrates from t_0 to t_f.
        %     X: the independent variables x (without t) as a cell array of 
        %        strings or as symbolic vector
        %        (optional: not needed if T is the only independent variable)
        %
        % Output:
        %     OBJ: a MODEL_IVP object with the passed configurations
        %
        
            % convert input
            f_sym = simplify(util.make_sym(f));
            p_sym = util.make_column_vector(util.make_sym(p));
            y_sym = util.make_sym(y);
            y0_sym = util.make_sym(y0);
            t_sym = util.make_sym(t);
            t_interval = util.make_column_vector(t_interval);
            if nargin >= 7
                x_sym = util.make_column_vector(util.make_sym(x));
            else
                x_sym = [];
            end
            if nargin < 8
                ode_solver_name = 'ode45';
            end
            this.ode_solver_name = ode_solver_name;
            
            % store input
            this.f_sym = f_sym;
            this.p_sym = p_sym;
            this.y_sym = y_sym;
            this.y0_sym = y0_sym;
            this.t_sym = t_sym;
            this.t_interval = t_interval;
            this.x_sym = x_sym;
            
            this.set_debug(0);
            
            % derivation function
            function dx_f = d(f, x)
                dx_f = simplify(jacobian(f, x));
            end
            
            n = length(p_sym);
            dp_y_sym = sym(['dp_' char(y_sym)], [1, n]);
            this.dp_y_sym = dp_y_sym;
            
            dp_f_sym = d(f_sym, y_sym) * dp_y_sym + d(f_sym, p_sym);
            this.dp_f_sym = dp_f_sym;
            
            dp_y0_sym = d(y0_sym, p_sym);
            this.dp_y0_sym = dp_y0_sym;
            
            dpdp_y_sym = sym(['dpdp_' char(y_sym)], [n, n]);
            this.dpdp_y_sym = dpdp_y_sym;
            
            dpdp_f_sym = d(dp_f_sym, dp_y_sym) * dpdp_y_sym + d(dp_f_sym, y_sym) * dp_y_sym + d(dp_f_sym, p_sym);
            dpdp_f_sym = simplify(dpdp_f_sym);
            this.dpdp_f_sym = dpdp_f_sym;
            
            dpdp_y0_sym = d(dp_y0_sym, p_sym);
            this.dpdp_y0_sym = dpdp_y0_sym;
            
            % add remove listeners
            addlistener(this, 'event_p_changed', @(src, evnt_data) remove_calculations('p'));
            addlistener(this, 'event_t_changed', @(src, evnt_data) remove_calculations('t'));
            addlistener(this, 'event_x_changed', @(src, evnt_data) remove_calculations('x'));
            
            function remove_calculations(changed)
                if strfind('p x',  changed)
                    this.y_p_x = [];
                    this.dp_y_p_x = [];
                    this.dpdp_y_p_x = [];
                end
                if strfind('p x t',  changed)
                    this.y_p_x_t = [];
                    this.dp_y_p_x_t = [];
                    this.dpdp_y_p_x_t = [];
                end
            end
        end
        
        
        function M = get_M(this, p, x_and_maybe_t)
        % GET_M returns the solution of the initial value problem with model
        % parameters P and the independent variables X and T.
        %
        % Example:
        %     M = MODEL_IVP_OBJECT.GET_M(P, X_AND_MAYBE_T)
        %
        % Input:
        %     P: the parameters
        %     X_AND_MAYBE_T: the current values for the independent variables X
        %         as a vector or the current values for the independent 
        %         variables X and T as one vector where the value of T is the
        %         last value
        %
        % Output:
        %     M: if T is passed: the solution of the initial value
        %        problem for the passed parameters P and the passed
        %        independent variables X and T
        %        if T is not passed: a handle to the solution function of
        %        the initial value problem for the passed parameters P and the
        %        passed independent variables X
        %
        
            % set input
            if nargin >= 2
                this.set_p(p);
            end
            if nargin >= 3
                this.set_x_and_maybe_t(x_and_maybe_t);
            end
            
            % solve initial value problem     
            if isempty(this.y_p_x)
                this.show_debug('calculating y');
                this.show_debug('   solving ODE 1 of 1');
            
                % get ode
                f_sym = this.get_f_sym();
                
                % substitude values for p
                p_sym = this.get_p_sym();
                p = this.get_p();
                f_p_sym= simplify(subs(f_sym, p_sym, p));
                
                % substitude values for x
                x_sym = this.get_x_sym();
                x = this.get_x();
                f_p_x_sym = simplify(subs(f_p_sym, x_sym, x));
                
                % prepare ode function
                t_sym = this.get_t_sym();
                y_sym = this.get_y_sym();
                f_p_x = @(t, y) (double(subs(subs(f_p_x_sym, [t_sym, y_sym], [t, y]))));
                
                % substitude values for p and x in y0
                y0_sym = this.get_y0_sym();
                y0_p_sym = subs(y0_sym, p_sym, p);
                y0_p_x_sym = subs(y0_p_sym, x_sym, x);
                y0_p_x = double(y0_p_x_sym);
                
                % solve ivp
                t_interval = this.get_t_interval();
                y_p_x = this.solve_ODE(f_p_x, t_interval, y0_p_x);
                
                % save result
                this.y_p_x = y_p_x;
                
                this.show_debug('   solving ODE 1 of 1 END');                    
                this.show_debug('calculating y END');
            else
                y_p_x = this.y_p_x;
            end
                        
            % apply t if passed
            if nargin >= 3 && length(x_and_maybe_t) > length(this.x_sym)
                if isempty(this.y_p_x_t)
                    t = this.get_t();
                    y_p_x_t = y_p_x(t);
                    this.y_p_x_t = y_p_x_t;
                else
                    y_p_x_t = this.y_p_x_t;
                end
                M = y_p_x_t;
            else
                M = y_p_x;
            end
            
        end
        
        function dp_M = get_dp_M(this, p, x_and_maybe_t)
        % GET_M returns the first derivative with respect to the parameters P of
        % the solution of the initial value problem with model parameters P and
        % the independent variables X and T.
        %
        % Example:
        %     DP_M = MODEL_IVP_OBJECT.GET_DP_M(P, X_AND_MAYBE_T)
        %
        % Input:
        %     P: the parameters
        %     X_AND_MAYBE_T: the current values for the independent variables X
        %         as a vector or the current values for the independent 
        %         variables X and T as one vector where the value of T is the
        %         last value
        %
        % Output:
        %     DP_M: if T is passed: the first derivative with respect to the
        %           parameters of the solution of the initial value
        %           problem for the passed parameters P and the passed
        %           independent variables X and T
        %           if X is not passed: a handle to the first derivative with
        %           respect to the parameters P of the solution function of
        %           the initial value problem for the passed parameters P and
        %           the passed independent variables X
        %
        
            % set input
            if nargin >= 2
                this.set_p(p);    
            end
            if nargin >= 3
                this.set_x_and_maybe_t(x_and_maybe_t);
            end
            
            % calculate output
            if isempty(this.dp_y_p_x)
                this.show_debug('calculating dp_y');
                
                % get ode
                dp_f_sym = this.get_dp_f_sym();
                
                % substitude values for p
                p_sym = this.get_p_sym();
                p = this.get_p();
                dp_f_p_sym= simplify(subs(dp_f_sym, p_sym, p));
                
                % substitude values for x
                x_sym = this.get_x_sym();
                x = this.get_x();
                dp_f_p_x_sym = simplify(subs(dp_f_p_sym, x_sym, x));
                
                % substitude values for p and x in dp_y0
                dp_y0_sym = this.get_dp_y0_sym();
                dp_y0_p_sym = subs(dp_y0_sym, p_sym, p);
                dp_y0_p_x_sym = subs(dp_y0_p_sym, x_sym, x);
                dp_y0_p_x = double(dp_y0_p_x_sym);

                % solve ivp
                y_p_x = this.get_M();

                % prepare to solve derivative ivps
                t_sym = this.get_t_sym();
                y_sym = this.get_y_sym();
                dp_y_sym = this.get_dp_y_sym();
                t_interval = this.get_t_interval();
                
                n = length(p);
                dp_y_p_x_cell = cell(n, 1);
                
                % solve derivative ivps in parallel
                parfor i=1:n
                    this.show_debug(['   solving ODE ' num2str(i) ' of ' num2str(n)]);
                    
                    dpi_f_p_x = @(t, dp_y) (double(subs(subs(dp_f_p_x_sym(i), [t_sym, y_sym, dp_y_sym(i)], [t, y_p_x(t), dp_y]))));
                    dp_y_p_x_cell{i} = this.solve_ODE(dpi_f_p_x, t_interval, dp_y0_p_x(i));
                    
                    this.show_debug(['   solving ODE ' num2str(i) ' of ' num2str(n) ' END']);
                end
                
                % make callable solution
                dp_y_p_x = @(t) eval_grad_cells(dp_y_p_x_cell, t);
       
                % save result
                this.dp_y_p_x = dp_y_p_x;
                
                this.show_debug('calculating dp_y END');
            else
                dp_y_p_x = this.dp_y_p_x;
            end
            
            function grad = eval_grad_cells(grad_cells, t)
                m = length(grad_cells);
                grad = zeros(m, 1);
                parfor i=1:m
                    grad(i) = grad_cells{i}(t);
                end
            end

            % apply t if passed
            if nargin >= 3 && length(x_and_maybe_t) > length(this.x_sym)
                if isempty(this.dp_y_p_x_t)
                    t = this.get_t();
                    dp_y_p_x_t = dp_y_p_x(t);
                    this.dp_y_p_x_t = dp_y_p_x_t;
                else
                    dp_y_p_x_t = this.dp_y_p_x_t;
                end
                dp_M = dp_y_p_x_t;
            else
                dp_M = dp_y_p_x;
            end
            
        end
        
        function dpdp_M = get_dpdp_M(this, p, x_and_maybe_t)
        % GET_M returns the second derivative with respect to the parameters P of
        % the solution of the initial value problem with model parameters P and
        % the independent variables X and T.
        %
        % Example:
        %     DPDP_M = MODEL_IVP_OBJECT.GET_DP_M(P, X_AND_MAYBE_T)
        %
        % Input:
        %     P: the parameters
        %     X_AND_MAYBE_T: the current values for the independent variables X
        %         as a vector or the current values for the independent 
        %         variables X and T as one vector where the value of T is the
        %         last value
        %
        % Output:
        %     DPDP_M: if T is passed: the second derivative with respect to the
        %           parameters of the solution of the initial value
        %           problem for the passed parameters P and the passed
        %           independent variables X and T
        %           if X is not passed: a handle to the second derivative with
        %           respect to the parameters P of the solution function of
        %           the initial value problem for the passed parameters P and
        %           the passed independent variables X
        %
        
            % set input
            if nargin >= 2
                this.set_p(p);    
            end
            if nargin >= 3
                this.set_x_and_maybe_t(x_and_maybe_t);
            end
            
            % calculate output
            if isempty(this.dpdp_y_p_x)
                this.show_debug('calculating dp_dp_y');
                
                % get ode
                dpdp_f_sym = this.get_dpdp_f_sym();
                
                % substitude values for p
                p_sym = this.get_p_sym();
                p = this.get_p();
                dpdp_f_p_sym= simplify(subs(dpdp_f_sym, p_sym, p));
                
                % substitude values for x
                x_sym = this.get_x_sym();
                x = this.get_x();
                dpdp_f_p_x_sym = simplify(subs(dpdp_f_p_sym, x_sym, x));
                
                % substitude values for p and x in dpdp_y0
                dpdp_y0_sym = this.get_dpdp_y0_sym();
                dpdp_y0_p_sym = subs(dpdp_y0_sym, p_sym, p);
                dpdp_y0_p_x_sym = subs(dpdp_y0_p_sym, x_sym, x);
                dpdp_y0_p_x = double(dpdp_y0_p_x_sym);

                % solve ivp
                y_p_x = this.get_M();
                dp_y_p_x = this.get_dp_M();

                % prepare to solve derivative ivps
                t_sym = this.get_t_sym();
                y_sym = this.get_y_sym();
                dp_y_sym = this.get_dp_y_sym();
                dpdp_y_sym = this.get_dpdp_y_sym();
                t_interval = this.get_t_interval();
                
                % prepare cell array and indices
                n = length(p);
                l = n*(n+1)/2;
                dpdp_y_p_x_cell = cell(l, 1);
                
                indices = zeros(l, 2);
                for k = 1:l
                    [i, j] = index_1d_to_2d(k);
                    indices(k, 1) = i;
                    indices(k, 2) = j;
                end
                
                % solve derivative derivative ivps in parallel
                parfor k = 1:l
                    i = indices(k, 1);
                    j = indices(k, 2);
                    
                    this.show_debug(['   solving ODE ' num2str(k) ' of ' num2str(l)]);  
                    
                    dpidpj_f_p_x = @(t, dpdp_y) (double(subs(subs(dpdp_f_p_x_sym(i, j), [t_sym, y_sym, dp_y_sym, dpdp_y_sym(i, j)], [t, y_p_x(t), dp_y_p_x(t).', dpdp_y]))));
                    dpdp_y_p_x_cell{k} = this.solve_ODE(dpidpj_f_p_x, t_interval, dpdp_y0_p_x(i, j));
                    
                    this.show_debug(['   solving ODE ' num2str(k) ' of ' num2str(l) ' END']);                    
                end
                
                % make callable solution
                dpdp_y_p_x = @(t) eval_hesse_cells(dpdp_y_p_x_cell, t);
       
                % save result
                this.dpdp_y_p_x = dpdp_y_p_x;
                
                this.show_debug('calculating dp_dp_y END');
            else
                dpdp_y_p_x = this.dpdp_y_p_x;
            end
            
            function [i, j] = index_1d_to_2d(k)
                i = ceil(((1+8*k)^(1/2) - 1) / 2);
                j = k - (i-1) * i / 2;
            end
            
            function hesse = eval_hesse_cells(hesse_cells, t)
                l = length(hesse_cells);
                n = index_1d_to_2d(l);
                hesse = zeros(n, n);

                for k=1:l
                    [i,j] = index_1d_to_2d(k);
                    hesse(i,j) = hesse_cells{k}(t);
                    hesse(j,i) = hesse(i,j);
                end
            end

            % apply t if passed
            if nargin >= 3 && length(x_and_maybe_t) > length(this.x_sym)
                if isempty(this.dpdp_y_p_x_t)
                    t = this.get_t();
                    dpdp_y_p_x_t = dpdp_y_p_x(t);
                    this.dpdp_y_p_x_t = dpdp_y_p_x_t;
                else
                    dpdp_y_p_x_t = this.dpdp_y_p_x_t;
                end
                dpdp_M = dpdp_y_p_x_t;
            else
                dpdp_M = dpdp_y_p_x;
            end

        end
        
        
        function set_debug(this, debug)
        % SET_DEBUG enables or disables the debug output.
        %
        % Example:
        %     MODEL_IVP_OBJECT.SET_DEBUG(DEBUG)
        %
        % Input:
        %     DEBUG: zero for disabling the debug output
        %            a positive integer for enabling the debug output
        %
        
            this.debug = debug;
        end
        
    end
    
    
    methods (Access = protected)
        
        function y_sol = solve_ODE(this, f, t_interval, y0, opt)
        % SOLVE_ODE solves an initial value problem.
        %
        % Example:
        %     Y_SOL = MODEL_IVP.SOLVE_ODE(F, T_interval, Y0)
        %
        % Input:
        %     F: a handle to the function for the derivative of Y with respect to X. F can
        %        depend on X and Y. In other words it holds (dy/dx)(x) = f(x,y).
        %     T_interval: a vector specifying the interval of integration, T_interval = [X0, Xf]. 
        %             The solver imposes the initial conditions at X0, and
        %             integrates from X0 to Xf.
        %     Y0: the initial value of Y at X0. In other words it holds
        %         y(x0) = y0.
        %
        % Output:
        %     Y_SOL: a handle to the solution function of the initial value problem.
        %
        
            % prepare default options if not passed
            if nargin < 5
                opt = odeset('RelTol', 10^-3, 'AbsTol', 10^-6);
                if this.use_debug(1)
                    opt.Stats = 'on';
                end
            end
            
            % chose solver
            ode_solver_name = this.ode_solver_name;
            ode_solver = str2func(ode_solver_name);
            
            % make tolerance not met warning catchable
            msg_id =  ['MATLAB:' ode_solver_name ':IntegrationTolNotMet'];
            warning('error', msg_id);
            
            % solve ode (and increase tolerance if unmetable)
            solved = false;
            while ~solved
                try
                    y_sol = ode_solver(f, t_interval, y0, opt);
                    solved = true;
                catch err
                   if (strcmp(err.identifier, msg_id))
                        opt.RelTol = opt.RelTol * 10;
                        opt.AbsTol = opt.AbsTol * 10^2;
                        warning(this.get_message_identifier('solve_ODE', 'integration_tolerances_relaxed'), 'Unable to meet integration tolerances without reducing the step size below the smallest value allowed. For this reason, relaxing integration tolerance.');
                   else
                      warning('on', msg_id);
                      rethrow(err);
                   end
                end
            end
            
            % make tolerance not met warning uncatchable
            warning('on', msg_id);
            
            % return result
            y_sol = @(t) deval(y_sol, t);
        end
        
        
        function set_p(this, p)
        % SET_P sets the model parameters P.
        %
        % Example:
        %     MODEL_IVP.SET_P(P)
        %
        % Input:
        %     P: the model parameters P
        %
        % see also GET_P
        %
        
            p = util.make_column_vector(p);
            p_sym = this.get_p_sym();
            
            if not(isequal(size(p), size(p_sym)))
                error(this.get_message_identifier('set_p', 'wrong_dimension'), 'The dimensions of p and p_sym must be the same.');
            end
            
            if not(isequal(this.p, p))
                this.p = p;
                notify(this, 'event_p_changed');
            end            
        end
        
        function p = get_p(this)
        % GET_P returns the model parameters P.
        %
        % Example:
        %     P = MODEL_IVP.GET_P()
        %
        % Output:
        %     P: the model parameters P
        %
        % see also SET_P
        %
        
            if isempty(this.p)
                error(this.get_message_identifier('get_p', 'p_not_set'), 'p is not set.');
            else
                p = this.p;
            end 
        end
         
        function p_sym = get_p_sym(this)
        % GET_P_SYM returns the symbolic variables representing the model
        % parameters P as a vector.
        %
        % Example:
        %     P_SYM = MODEL_IVP.GET_P_SYM()
        %
        % Output:
        %     P_SYM: the symbolic variables representing the model parameters P
        %
        
            p_sym = this.p_sym;
        end
        
        
        function set_x_and_maybe_t(this, x_and_maybe_t)
        % SET_X_AND_MAYBE_T sets the current values for the independent variables
        % X and T.
        %
        % Example:
        %     MODEL_IVP.SET_X_AND_MAYBE_T(X_AND_MAYBE_T)
        %
        % Input:
        %     X_AND_MAYBE_T: the current values for the independent variables X
        %         as a vector or the current values for the independent 
        %         variables X and T as one vector where the value of T is the
        %         last value
        %
        % see also GET_X, GET_T
        %
            x_and_maybe_t = util.make_column_vector(x_and_maybe_t);
            x_and_maybe_t_len = length(x_and_maybe_t);

            x_sym = this.get_x_sym();
            x_sym_len = length(x_sym);

            % check length
            if x_and_maybe_t_len < x_sym_len || x_and_maybe_t_len > x_sym_len + 1
                error(this.get_message_identifier('set_x_and_maybe_t', 'wrong_dimension'), 'The length of set_x_and_maybe_t must be the length of x_sym or the length of x_sym + 1.');
            end
            
            % set x
            if x_sym_len > 0
                if x_and_maybe_t_len == x_sym_len
                    x = x_and_maybe_t;
                else
                    x = x_and_maybe_t(1 : x_sym_len);
                end
                this.set_x(x);
            end
            
            % set t
            if x_and_maybe_t_len == x_sym_len + 1
                t = x_and_maybe_t(x_sym_len + 1);
                this.set_t(t);
            end
        end
        
        
        function set_x(this, x)
        % SET_X sets the current values for the independent variables X.
        %
        % Example:
        %     MODEL_IVP.SET_X(X)
        %
        % Input:
        %     X: the current values for the independent variables X
        %
        % see also GET_X
        %
            x = util.make_column_vector(x);
            x_sym = this.get_x_sym();
            
            if not(isequal(size(x), size(x_sym)))
                error(this.get_message_identifier('set_x', 'wrong_dimension'), 'The dimensions of x and x_sym must be the same.');
            end
            
            if not(isequal(this.x, x))
                this.x = x;
                notify(this, 'event_x_changed');
            end            
        end
        
        function x = get_x(this)
        % GET_X returns the current values for the independent variables X.
        %
        % Example:
        %     X = MODEL_IVP.GET_X()
        %
        % Output:
        %     X: the current values for the independent variables X
        %
        % see also SET_X
        %
        
            if ~ isempty(this.get_x_sym) && isempty(this.x)
                error(this.get_message_identifier('get_x', 'not_set'), 'x is not set.');
            else
                x = this.x;
            end 
        end
        
        function x_sym = get_x_sym(this)
        % GET_X_SYM returns the symbolic variables representing the independent
        % model variables X as a vector.
        %
        % Example:
        %     X_SYM = MODEL_IVP.GET_X_SYM()
        %
        % Output:
        %     X_SYM: the symbolic variables representing the independent model
        %            variables X
        %
        
            x_sym = this.x_sym;
        end  
        
        
        function set_t(this, t)
        % SET_T sets the current values for the independent variable T.
        %
        % Example:
        %     MODEL_IVP.SET_T(T)
        %
        % Input:
        %     T: the current values for the independent variable T
        %
        % see also GET_T
        %

            if not(isscalar(t))
                error(this.get_message_identifier('set_t', 'wrong_dimension'), 'T has to be a scalar.');
            end
            
            if not(isequal(this.t, t))
                this.t = t;
                notify(this, 'event_t_changed');
            end            
        end
        
        function t = get_t(this)
        % GET_T returns the current values for the independent variable T.
        %
        % Example:
        %     T = MODEL_IVP.GET_T()
        %
        % Output:
        %     T: the current values for the independent variable T
        %
        % see also SET_T
        %
        
            if isempty(this.t)
                error(this.get_message_identifier('get_t', 'not_set'), 'T is not set.');
            else
                t = this.t;
            end 
        end
        
        function t_sym = get_t_sym(this)
        % GET_T_SYM returns the symbolic variables representing the independent
        % model variables T of the ordinary differential 
        %
        % Example:
        %     T_SYM = MODEL_IVP.GET_T_SYM()
        %
        % Output:
        %     T_SYM: the symbolic variables representing the independent model
        %            variables T
        %
        
            t_sym = this.t_sym;
        end  
        
        
        function t_interval = get_t_interval(this)
        % GET_T_interval returns a vector specifying the interval of integration.
        %
        % Example:
        %     T_interval = MODEL_IVP.GET_T_interval()
        %
        % Output:
        %     T_interval: a vector specifying the interval of integration, T_interval = [X0, Xf]. 
        %             The solver imposes the initial conditions at X0, and
        %             integrates from X0 to Xf.
        %
        
            t_interval = this.t_interval;
        end
        
                
        function f_sym = get_f_sym(this)
        % GET_F_SYM returns the symbolic formula for the derivative of Y with respect to X.
        %
        % Example:
        %     F_SYM = MODEL_IVP.GET_F_SYM()
        %
        % Output:
        %     F_SYM: the symbolic formula for the derivative of Y with respect
        %            to X. F can depend on X, Y and P. In other words it
        %            holds (dy/dx)(x,p) = f(x,y,p). 
        %
        
            f_sym = this.f_sym;
        end
        
        function y_sym = get_y_sym(this)
        % GET_Y_SYM returns the symbolic variable representing the dependent
        % model variable Y.
        %
        % Example:
        %     Y_SYM = MODEL_IVP.GET_Y_SYM()
        %
        % Output:
        %     Y_SYM: the dependent symbolic variable Y
        %
        
            y_sym = this.y_sym;
        end
        
        function y0_sym = get_y0_sym(this)
        % GET_Y0_SYM returns the initial value of Y at X0 as symbolic formula.
        %
        % Example:
        %     Y0_SYM = MODEL_IVP.GET_Y0_SYM()
        %
        % Output:
        %     Y0_SYM: the initial value of Y at X0 as symbolic formula. Y0_SYM 
        %             can depend on the parameters P. In other words it holds
        %             y(x0,p) = y0(p).
        %
        
            y0_sym = this.y0_sym;
        end        
        
        
        function dp_f_sym = get_dp_f_sym(this)
        % GET_DP_F_SYM returns the symbolic formula for the derivative of Y with respect to X and P.
        %
        % Example:
        %     DP_F_SYM = MODEL_IVP.GET_DP_F_SYM()
        %
        % Output:
        %     DP_F_SYM: the symbolic formula for the derivative of Y with respect
        %            to X and P. DP_F can depend on X, Y and P. In other words 
        %            it holds (dy^2/(dx dp))(x,p) = dp_f(x,y,p). 
        %
        
            dp_f_sym = this.dp_f_sym;
        end
        
        function dp_y_sym = get_dp_y_sym(this)
        % GET_DP_Y_SYM returns the derivative of Y with respect to the parameters P as symbolic variable.
        %
        % Example:
        %     DP_Y_SYM = MODEL_IVP.GET_DP_Y_SYM()
        %
        % Output:
        %     DP_Y_SYM: the derivative of Y with respect to the parameters P
        %               as symbolic variable
        %
        
            dp_y_sym = this.dp_y_sym;
        end
        
        function dp_y0_sym = get_dp_y0_sym(this)
        % GET_DP_Y0_SYM returns the initial value of the derivative of Y with respect to P at X0 as symbolic formula.
        %
        % Example:
        %     DP_Y0_SYM = MODEL_IVP.GET_DP_Y0_SYM()
        %
        % Output:
        %     DP_Y0_SYM: the initial value of the derivative of Y with respect
        %                to P at X0 as symbolic formula DP_Y0_SYM can
        %                depend on the parameters P. In other words it
        %                holds (dy/dp)(x0,p) = dp_y0(p).
        %
        
            dp_y0_sym = this.dp_y0_sym;
        end
        
        
        function dpdp_f_sym = get_dpdp_f_sym(this)
        % GET_DPDP_F_SYM returns the symbolic formula for the derivative of Y with respect to once X and twice P.
        %
        % Example:
        %     DPDP_F_SYM = MODEL_IVP.GET_DPDP_F_SYM()
        %
        % Output:
        %     DPDP_F_SYM: the symbolic formula for the derivative of Y with
        %                 respect to once X and twice P. DP_F can depend on
        %                 X, Y and P. In other words it holds 
        %                 (dy^3/(dx dp^2))(x,p) = dpdp_f(x,y,p). 
        %
        
            dpdp_f_sym = this.dpdp_f_sym;
        end
        
        function dpdp_y_sym = get_dpdp_y_sym(this)
        % GET_DPDP_Y_SYM returns the second derivative of Y with respect to the parameters P as symbolic variable.
        %
        % Example:
        %     DPDP_Y_SYM = MODEL_IVP.GET_DPDP_Y_SYM()
        %
        % Output:
        %     DPDP_Y_SYM: the second derivative of Y with respect to the parameters P
        %                 as symbolic variable
        %
        
            dpdp_y_sym = this.dpdp_y_sym;
        end
        
        function dpdp_y0_sym = get_dpdp_y0_sym(this)
        % GET_DPDP_Y0_SYM returns the initial value of the second derivative of Y with respect to P at X0 as symbolic formula.
        %
        % Example:
        %     DPDP_Y0_SYM = MODEL_IVP.GET_DPDP_Y0_SYM()
        %
        % Output:
        %     DPDP_Y0_SYM: the initial value of the second derivative of Y
        %                  with respect to P at X0 as symbolic formula
        %                  DP_Y0_SYM can depend on the parameters P. In
        %                  other words it holds (dy^2/dp^2)(x0,p) = dpdp_y0(p).
        %
        
            dpdp_y0_sym = this.dpdp_y0_sym;
        end
        
        
        function use_debug = use_debug(this, debug)
        % USE_DEBUG returns whether debug informations will be output or not.
        %
        % Example:
        %     MODEL_IVP.USE_DEBUG(DEBUG)
        %
        % Input:
        %     DEBUG: the debug-level which will be checked
        %            (optional, default: 1)
        %
        % Output:
        %     USE_DEBUG: return whether the current debug-level is greater or equal the passed debug-level
        %
        % see also SET_DEBUG
        %
        
            if nargin <= 1
                debug = 1;
            end
            use_debug = this.debug >= debug;
        end
        
                
        function show_debug(this, message)
        % SHOW_DEBUG outputs a debug information with the passed message.
        %
        % Example:
        %     MODEL_IVP.SHOW_DEBUG(MESSAGE)
        %
        % Input:
        %     MESSAGE: the message for the debug information
        %
            if this.use_debug(2)
                display(['MODEL_IVP: ' message]);
            end
        end
        
    end
    
    methods (Access = protected, Static)
        
        function s = get_message_identifier(method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in methods of these object.
        %
        % Example:
        %     ID = MODEL_IVP.GET_MESSAGE_IDENTIFIER(METHOD, MNEMONIC)
        %
        % Input:
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
            s = util.get_message_identifier('model_ivp', method, mnemonic);
        end       
        
    end
        
end
    


%#ok<*PROP>
