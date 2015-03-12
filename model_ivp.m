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
    Copyright (C) 2010-2015 Joscha Reimer jor@informatik.uni-kiel.de

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
        x_sym;
        p_sym;
        
        x_span;
        p;
        x;
        
        y_sym;
        f_sym;
        y0_sym;
        y;
        yx;
        
        dp_y_sym;
        dp_f_sym;
        dp_y0_sym;    
        dp_y;
        dp_yx;
        
        dpdp_y_sym;
        dpdp_f_sym;
        dpdp_y0_sym;
        dpdp_y;
        dpdp_yx;
        
        debug;
    end
    
    events (ListenAccess = protected, NotifyAccess = protected)
        event_p_changed;    %is triggered when p was changed
        event_x_changed;    %is triggered when t was changed
    end
    
    methods (Access = public)
        
        function this = model_ivp(f, p, y, y0, x, x_span)
        % MODEL_IVP creates an MODEL_IVP object.
        %
        % Example:
        %     OBJ = MODEL_IVP(F, P, Y, Y0, X, X_SPAN)
        %
        % Input:
        %     F: the formula f of the derivative of y with respect to x as
        %        string or symbolic formula. f can depend on x, y and p.
        %        It holds (dy/dx)(x,p) = f(x,y,p). 
        %     P: the variables of the parameters p as a cell array of
        %        strings or as symbolic vector
        %     Y: the independent variable y as string or symbolic variable
        %     Y0: the value of y at x0 as string or symbolic variable.
        %         Y0_SYM can depend on the parameters p. It holds y(x0,p) = y0(p).
        %     X: the dependent variable x as string or symbolic variable
        %     X_SPAN: a vector specifying the interval of integration, X_SPAN = [x0, xf]. 
        %             The solver imposes the initial conditions at x0, and
        %             integrates from x0 to xf.
        %
        % Output:
        %     OBJ: a MODEL_IVP object with the passed configurations
        %
        
            x_sym = util.make_sym(x);
            y_sym = util.make_sym(y);
            p_sym = util.make_column_vector(util.make_sym(p));
            f_sym = util.make_sym(f);
            y0_sym = util.make_sym(y0);
            x_span = util.make_column_vector(x_span);
            this.x_sym = x_sym;
            this.y_sym = y_sym;
            this.p_sym = p_sym;
            this.f_sym = f_sym;
            this.y0_sym = y0_sym;
            this.x_span = x_span;
            
            
            % prepare dt_dp_y function
            for i=length(p_sym):-1:1
                dp_y_sym(i) = sym(['d' char(p_sym(i)) '_y']);
            end
            this.dp_y_sym = dp_y_sym;
            
            dp_f_sym = d(f_sym, y_sym) * dp_y_sym + d(f_sym, p_sym);
            this.dp_f_sym = dp_f_sym;
            
            dp_y0_sym = d(y0_sym, p_sym);
            this.dp_y0_sym = dp_y0_sym;
            
            
            % prepare dt_dp_dp_y function
            for i=length(p_sym):-1:1
                for j=length(p_sym):-1:1
                    dpdp_y_sym(i,j) = sym(['d' char(p_sym(i)) 'd' char(p_sym(j)) '_y']);
                end
            end
            this.dpdp_y_sym = dpdp_y_sym;
            
            dpdp_f_sym = d(dp_f_sym, dp_y_sym) * dpdp_y_sym + d(dp_f_sym, y_sym) * dp_y_sym + d(dp_f_sym, p_sym);
            dpdp_f_sym = simplify(dpdp_f_sym);
            this.dpdp_f_sym = dpdp_f_sym;
            
            dpdp_y0_sym = d(dp_y0_sym, p_sym);
            this.dpdp_y0_sym = dpdp_y0_sym;
            
            this.set_debug(0);
            
            function dx_f = d(f, x)
                dx_f = simplify(jacobian(f, x));
            end
            
            
            addlistener(this, 'event_p_changed', @(src, evnt_data) remove_calculations('p'));
            addlistener(this, 'event_x_changed', @(src, evnt_data) remove_calculations('t'));
            
            function remove_calculations(changed)
                if strfind('p',  changed)
                    this.y = [];
                    this.dp_y = [];
                    this.dpdp_y = [];
                end
                if strfind('p t',  changed)
                    this.yx = [];
                    this.dp_yx = [];
                    this.dpdp_yx = [];
                end
            end
        end
        
        
        function M = get_M(this, p, x)
        % GET_M returns the solution of the initial value problem with parameter P.
        %
        % Example:
        %     M = MODEL_IVP_OBJECT.GET_M(P, X)
        %
        % Input:
        %     P: the parameters
        %     X: the dependent variable (optional)
        %
        % Output:
        %     M: if X is passed: the solution of the initial value
        %        problem for the passed parameters P and the passed
        %        dependent variable X
        %        if X is not passed: a handle to the solution function of
        %        the initial value problem for the passed parameters P
        %
        
            if nargin >= 2
                this.set_p(p);    
            end
            if nargin >= 3
                this.set_x(x);
            end
                        
            if isempty(this.y)
                this.show_debug('calculating y');                    
                this.show_debug('   solving ODE 1 of 1');
            
                f_sym = this.get_f_sym();
                p_sym = this.get_p_sym();
                p = this.get_p();
                x_sym = this.get_x_sym();            
                y_sym = this.get_y_sym();
                y0_sym = this.get_y0_sym();
                x_span = this.get_x_span(); 
                
                f_tmp = simplify(subs(f_sym, p_sym, p));
                f = @(x, y) (double(subs(subs(f_tmp, [x_sym, y_sym], [x, y]))));
                               
                y0 = double(subs(y0_sym, p_sym, p));
                
                y = this.solve_ODE(f, x_span, y0);
                
                this.y = y;
                
                this.show_debug('   solving ODE 1 of 1 END');                    
                this.show_debug('calculating y END');
            else
                y = this.y;
            end
                        
            
            if nargin >= 3
                x = this.get_x();
                
                if isempty(this.yx)
                    yx = y(x);
                    this.yx = yx;
                else
                    yx = this.yx;
                end
                M = yx;
            else
                M = y;
            end
            
        end
        
        function dp_M = get_dp_M(this, p, x)
        % GET_DP_M returns the first derivative with respect to the parameters P of the solution of the initial value problem.
        %
        % Example:
        %     DP_M = MODEL_IVP_OBJECT.GET_DP_M(P, X)
        %
        % Input:
        %     P: the parameters
        %     X: the dependent variable (optional)
        %
        % Output:
        %     DP_M: if X is passed: the first derivative with respect to the
        %           parameters of the solution of the initial value
        %           problem for the passed parameters P and the passed
        %           dependent variable X
        %           if X is not passed: a handle to the first derivative with
        %           respect to the parameters P of the solution function of
        %           the initial value problem for the passed parameters P
        %
        
            if nargin >= 2
                this.set_p(p);    
            end
            if nargin >= 3
                this.set_x(x);
            end
                       
            if isempty(this.dp_y)
                this.show_debug('calculating dp_y');
                
                dp_f_sym = this.get_dp_f_sym();
                dp_y0_sym = this.get_dp_y0_sym();
                p_sym = this.get_p_sym();
                p = this.get_p();
                
                dp_f_tmp = simplify(subs(dp_f_sym, p_sym, p));
                dp_y0 = double(subs(dp_y0_sym, p_sym, p));
                
                x_sym = this.get_x_sym();            
                y_sym = this.get_y_sym();
                dp_y_sym = this.get_dp_y_sym();           
                x_span = this.get_x_span();
                y = this.get_M(p);
                
                n = length(p);
                dp_y_cell = cell(n, 1);
                
                parfor i=1:n
                    this.show_debug(['   solving ODE ' num2str(i) ' of ' num2str(n)]);
                    
                    dpi_f = @(x, dp_y) (double(subs(subs(dp_f_tmp(i), [x_sym, y_sym, dp_y_sym(i)], [x, y(x), dp_y]))));
                    dp_y_cell{i} = this.solve_ODE(dpi_f, x_span, dp_y0(i));
                    
                    this.show_debug(['   solving ODE ' num2str(i) ' of ' num2str(n) ' END']);
                end
                
                dp_y = @(t) eval_grad_cells(dp_y_cell, t);
                this.dp_y = dp_y;                
                
                this.show_debug('calculating dp_y END');
            else
                dp_y = this.dp_y;
            end
            
            
            function grad = eval_grad_cells(grad_cells, x)
                m = length(grad_cells);
                grad = zeros(m, 1);
                for i=1:m
                    grad(i) = grad_cells{i}(x);
                end
                %{
                    n = length(t);
                    m = length(grad_cells);
                    grad = zeros(n, m);
                    for i=1:m
                        grad(:, i) = grad_cells{i}(t).';
                    end
                %}
            end
            
            
            if nargin >= 3
                x = this.get_x();
                
                if isempty(this.dp_yx)
                    dp_yx = dp_y(x);
                    this.dp_yx = dp_yx;
                else
                    dp_yx = this.dp_yx;
                end
                dp_M = dp_yx;
            else
                dp_M = dp_y;
            end
            
        end
        
        function dpdp_M = get_dpdp_M(this, p, x)
        % GET_DPDP_M returns the second derivative with respect to the parameters P of the solution of the initial value problem.
        %
        % Example:
        %     DPDP_M = MODEL_IVP_OBJECT.GET_DPDP_M(P, X)
        %
        % Input:
        %     P: the parameters
        %     X: the dependent variable (optional)
        %
        % Output:
        %     DPDP_M: if X is passed: the second derivative with respect to
        %             the parameters of the solution of the initial value
        %             problem for the passed parameters P and the passed
        %             dependent variable X
        %             if X is not passed: a handle to the second derivative
        %             with respect to the parameters P of the solution 
        %             function of the initial value problem for the passed
        %             parameters P
        %
        
            if nargin >= 2
                this.set_p(p);    
            end
            if nargin >= 3
                this.set_x(x);
            end
            
            if isempty(this.dpdp_y)
                this.show_debug('calculating dp_dp_y');                
                
                dpdp_f_sym = this.get_dpdp_f_sym();
                dpdp_y0_sym = this.get_dpdp_y0_sym();
                p_sym = this.get_p_sym();
                p = this.get_p();
                                
                dpdp_f_tmp = simplify(subs(dpdp_f_sym, p_sym, p));
                dpdp_y0 = double(subs(dpdp_y0_sym, p_sym, p));
                
                x_sym = this.get_x_sym();            
                y_sym = this.get_y_sym();   
                dp_y_sym = this.get_dp_y_sym();   
                dpdp_y_sym = this.get_dpdp_y_sym();   
                x_span = this.get_x_span();
                
                y = this.get_M(p);
                dp_y = this.get_dp_M(p);
                                
                n = length(p);
                l = n*(n+1)/2;
                dpdp_y_cell = cell(l, 1);
                
                indices = zeros(l, 2);
                for k = 1:l
                    [i, j] = index_1d_to_2d(k);
                    indices(k, 1) = i;
                    indices(k, 2) = j;
                end
                
                parfor k = 1:l
                    i = indices(k, 1);
                    j = indices(k, 2);
                    
                    this.show_debug(['   solving ODE ' num2str(k) ' of ' num2str(l)]);  
                    
                    dpidpj_f = @(x, dpdp_y) (double(subs(subs(dpdp_f_tmp(i, j), [x_sym, y_sym, dp_y_sym, dpdp_y_sym(i, j)], [x, y(x), dp_y(x).', dpdp_y]))));
                    dpdp_y_cell{k} = this.solve_ODE(dpidpj_f, x_span, dpdp_y0(i, j));
                    
                    this.show_debug(['   solving ODE ' num2str(k) ' of ' num2str(l) ' END']);                    
                end
                
                dpdp_y = @(x) eval_hesse_cells(dpdp_y_cell, x);
                this.dpdp_y = dpdp_y;
                
                this.show_debug('calculating dp_dp_y END');
            else
                dpdp_y = this.dpdp_y;
            end
            
            
            function [i, j] = index_1d_to_2d(k)
                i = ceil(((1+8*k)^(1/2) - 1) / 2);
                j = k - (i-1) * i / 2;
            end
            
            function hesse = eval_hesse_cells(hesse_cells, x)
                l = length(hesse_cells);
                n = index_1d_to_2d(l);
                hesse = zeros(n, n);

                for k=1:l
                    [i,j] = index_1d_to_2d(k);
                    hesse(i,j) = hesse_cells{k}(x);
                    hesse(j,i) = hesse(i,j);
                end
            end
            
            
            if nargin >= 3
                x = this.get_x();
                
                if isempty(this.dpdp_yx)
                    dpdp_yx = dpdp_y(x);
                    this.dpdp_yx = dpdp_yx;
                else
                    dpdp_yx = this.dpdp_yx;
                end
                dpdp_M = dpdp_yx;
            else
                dpdp_M = dpdp_y;
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
        
        function y_sol = solve_ODE(this, f, x_span, y0, opt)
        % SOLVE_ODE solves an initial value problem.
        %
        % Example:
        %     Y_SOL = MODEL_IVP.SOLVE_ODE(F, X_SPAN, Y0)
        %
        % Input:
        %     F: a handle to the function for the derivative of Y with respect to X. F can
        %        depend on X and Y. In other words it holds (dy/dx)(x) = f(x,y).
        %     X_SPAN: a vector specifying the interval of integration, X_SPAN = [X0, Xf]. 
        %             The solver imposes the initial conditions at X0, and
        %             integrates from X0 to Xf.
        %     Y0: the initial value of Y at X0. In other words it holds
        %         y(x0) = y0.
        %
        % Output:
        %     Y_SOL: a handle to the solution function of the initial value problem.
        %
        
            if nargin < 5
                opt = odeset('RelTol', 10^-3, 'AbsTol', 10^-6);
                if this.use_debug(1)
                    opt.Stats = 'on';
                end
            end
            
            msg_id =  'MATLAB:ode23s:IntegrationTolNotMet';
            warning('error', msg_id);
            solved = false;
            while ~solved
                try
                    y_sol = ode23s(f, x_span, y0, opt);
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
            warning('on', msg_id);
            y_sol = @(t) deval(y_sol, t);
        end
        
        
        function set_p(this, p)
        % SET_P sets the parameters P.
        %
        % Example:
        %     MODEL_IVP.SET_P(P)
        %
        % Input:
        %     P: the parameters P
        %
        % see also GET_P
        %
        
            p = util.make_column_vector(p);
            p_sym = this.get_p_sym();
            
            if not(isequal(size(p), size(p_sym)))
                error(this.get_message_identifier('set_p', 'wrong_dimension'), 'The dimensions of p and p_sym don''t match.');
            end
            
            if not(isequal(this.p, p))
                this.p = p;
                notify(this, 'event_p_changed');
            end            
        end
        
        function p = get_p(this)
        % GET_P returns the parameters P.
        %
        % Example:
        %     P = MODEL_IVP.GET_P()
        %
        % Output:
        %     P: the parameters P
        %
        % see also SET_P
        %
        
            if isempty(this.p)
                error(this.get_message_identifier('get_p', 'p_not_set'), 'p is not set.');
            else
                p = this.p;
            end 
        end
        
        function set_x(this, t)
        % SET_X sets the dependent variable X.
        %
        % Example:
        %     MODEL_IVP.SET_X(X)
        %
        % Input:
        %     X: the dependent variable X
        %
        % see also GET_X
        %
        
            t = util.make_column_vector(t);
            if not(isequal(this.x, t))
                this.x = t;
                notify(this, 'event_x_changed');
            end            
        end
        
        function t = get_x(this)
        % GET_X returns the dependent variable X.
        %
        % Example:
        %     X = MODEL_IVP.GET_X()
        %
        % Output:
        %     X: the dependent variable X
        %
        % see also SET_X
        %
        
            if isempty(this.x)
                error(this.get_message_identifier('get_x', 'x_not_set'), 'x is not set.');
            else
                t = this.x;
            end 
        end
         
        
        function x_sym = get_x_sym(this)
        % GET_X_SYM returns the dependent symbolic variable X.
        %
        % Example:
        %     X_SYM = MODEL_IVP.GET_X_SYM()
        %
        % Output:
        %     X_SYM: the dependent symbolic variable X
        %
        
            x_sym = this.x_sym;
        end        
                
        function y_sym = get_y_sym(this)
        % GET_Y_SYM returns the independent symbolic variable Y.
        %
        % Example:
        %     Y_SYM = MODEL_IVP.GET_Y_SYM()
        %
        % Output:
        %     Y_SYM: the independent symbolic variable Y
        %
        
            y_sym = this.y_sym;
        end
        
        function p_sym = get_p_sym(this)
        % GET_P_SYM returns the parameters P as symbolic vector.
        %
        % Example:
        %     P_SYM = MODEL_IVP.GET_P_SYM()
        %
        % Output:
        %     P_SYM: the parameters P as symbolic vector
        %
        
            p_sym = this.p_sym;
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
        
        function x_span = get_x_span(this)
        % GET_X_SPAN returns a vector specifying the interval of integration.
        %
        % Example:
        %     X_SPAN = MODEL_IVP.GET_X_SPAN()
        %
        % Output:
        %     X_SPAN: a vector specifying the interval of integration, X_SPAN = [X0, Xf]. 
        %             The solver imposes the initial conditions at X0, and
        %             integrates from X0 to Xf.
        %
        
            x_span = this.x_span;
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
