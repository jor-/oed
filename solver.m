classdef solver < handle
% SOLVER allows to calculate and optimize the quality of experimental designs and allows to calculate a parameter estimation resulting from accomplished measurements.
%
% SOLVER Methods:
%   GET_OPTIMAL_PARAMETERS - returns a parameter estimation resulting 
%       from the weighted least square method applied to the accomplished
%       measurements.
%   GET_OPTIMAL_MEASUREMENTS - returns the optimal subset of measurements
%       of the selectable measurements.
%   GET_OPTIMAL_WEIGHTS - returns the optimal weights of the measurements.
%   GET_QUALITY - returns the quality resulting from the passed weights
%       of the measurements.
%   SET_MODEL - sets the model which will be used for the computations. 
%   SET_INITIAL_PARAMETER_ESTIMATION - sets the initial estimation of the
%       model parameter. 
%   SET_SELECTABLE_MEASUREMENTS - sets the selectable measurements and the
%       associated variances.
%   SET_ACCOMPLISHED_MEASUREMENTS - sets the accomplished measurements
%       including the measurements and the associated variances. 
%   SET_OPTION - sets the option for this SOLVER object.
%
%
%   To calculate the optimal measurements use the GET_OPTIMAL_WEIGHTS 
%   method. Use the GET_QUALITY method to calculate the quality of a set of
%   measurements. You can calculate an parameter estimation from the 
%   accomplished measurements with the GET_OPTIMAL_PARAMETERS method.
%   
%   You have to provide a model that describes the observed process and an
%   estimation of the parameters of the model. Further you have to provide
%   selectable measurements with the associated variance of the
%   measurement error and if available accomplished measurements. You can
%   provide these things via the corresponding SET methods.

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

    properties (Access = private)
        model;
        p0;
        p;
        p_lb;
        p_ub;
        
        t_var;
        v_var;
        
        t_fix;
        v_fix;
        eta;
        
        w_var;
        
        options;
        
        
        dp_Mt;
        dpdp_Mt;
        
        C;
        dp_C; 
        dw_C;
        
        H;
        dw_H;
        
        phi;
        dp_phi;
        dw_phi;
        dwdp_phi;
        
        C0;
        
        phiR;
        dw_phiR;
    end
    
    events (ListenAccess = protected, NotifyAccess = protected)
        event_model_changed;        %is triggered when the model MODEL was changed
        event_p0_changed;           %is triggered when the initial parameter estimation P0 was changed
        event_p_changed;            %is triggered when the parameter estimation P was changed
        
        event_t_changed;            %is triggered when the measurements T were changed
        event_v_changed;            %is triggered when the variances V were changed
        event_eta_changed;          %is triggered when the measurement results were changed
        event_w_changed;            %is triggered when the weights W were changed
        
        event_alpha_changed;        %is triggered when the alpha was changed
        event_criterion_changed;    %is triggered when the quality criterion was changed
        %event_solver_changed;       %is triggered when the solver was changed
        event_estimation_method_changed;   %is triggered when the estimation method of the real parameter was changed
    end
    
    
    methods (Access = public)
        
        % ******************* CONSTRUCTOR ******************* %
        function this = solver(model, p0, t_var, v_var)
        % SOLVER creates a SOLVER object.
        %
        % Example:
        %     OBJ = SOLVER(MODEL, P, T_VAR, V_VAR)
        %
        % Input:
        %     MODEL: the model represented as an object whose class
        %            implements the MODEL interface
        %            (optional, can later be set by SET_MODEL)
        %     p0: the inital estimation of the model parameters
        %         (optional, can later be set by SET_INITIAL_PARAMETER_ESTIMATION)
        %     T_VAR: the measurements for which the quality will be
        %            calculated or optimized
        %            (optional, can later be set by SET_ACCOMPLISHED_MEASUREMENTS)
        %            format: a n x m matrix where n the dimension of a
        %                    measurement and m is the number of
        %                    measurements
        %     V_VAR: the variances of the measurement errors associated with 
        %            these measurements
        %            (optional, can later be set by SET_ACCOMPLISHED_MEASUREMENTS)
        %            format: a vector of length m where m is the number of 
        %                    measurements
        %
        % Output:
        %     OBJ: a SOLVER object with the passed configurations
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION, SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if nargin >= 1
                this.set_model(model);
            end
            if nargin >= 2
                this.set_initial_parameter_estimation(p0);
            end
            if nargin == 3
                warning(this.get_message_identifier('solver', 'v_var_missing'), 't_var is set but v_var is missing therefor t_var will be ignored');
            end
            if nargin >= 4
                this.set_selectable_measurements(t_var, v_var);
            end
            
            this.options = solver_options();
            
            
            model_changed = 1;
            p0_changed = 2;
            p_changed = 3;
            t_changed = 4;
            v_changed = 5;
            eta_changed = 6;
            w_changed = 7;
            criterion_changed = 8;
            alpha_changed = 9;
            estimation_changed = 9;
            
            addlistener(this, 'event_model_changed', @(src, evnt_data) remove_calculations(model_changed));
            addlistener(this, 'event_p0_changed', @(src, evnt_data) remove_calculations(p0_changed));
            addlistener(this, 'event_p_changed', @(src, evnt_data) remove_calculations(p_changed));
            
            addlistener(this, 'event_t_changed', @(src, evnt_data) remove_calculations(t_changed));
            addlistener(this, 'event_v_changed', @(src, evnt_data) remove_calculations(v_changed));
            addlistener(this, 'event_eta_changed', @(src, evnt_data) remove_calculations(eta_changed));
            addlistener(this, 'event_w_changed', @(src, evnt_data) remove_calculations(w_changed));
            
            addlistener(this, 'event_criterion_changed', @(src, evnt_data) remove_calculations(criterion_changed));
            addlistener(this, 'event_alpha_changed', @(src, evnt_data) remove_calculations(alpha_changed));
            addlistener(this, 'event_estimation_method_changed', @(src, evnt_data) remove_calculations(estimation_changed));
            %addlistener(this, 'event_solver_changed', @(src, evnt_data) remove_calculations(solver_changed));
            
            function remove_calculations(changed)
                if changed <= t_changed     % m, p0, p, t
                    this.dp_Mt = [];
                    this.dpdp_Mt = [];
                end
                if changed <= v_changed     % m, p0, p, t, v
                    this.C0 = [];
                end
                if changed <= eta_changed   % m, p0, p, t, v, eta
                    this.p = [];
                end
                if changed <= w_changed     % m, p0, p, t, v, eta, w
                    this.C = [];
                    this.dp_C = [];
                    this.dw_C = [];
                    this.H = [];
                    this.dw_H = [];
                end
                if changed <= criterion_changed   % m, p0, p, t, v, eta, w, c
                    this.phi = [];
                    this.dp_phi = [];
                    this.dw_phi = [];
                    this.dwdp_phi = [];
                end
                if changed <= alpha_changed       % m, p0, p, t, v, eta, w, c, alpha, estimation
                    this.phiR = [];
                    this.dw_phiR = [];
                end  
            end
            
        end
                
        
        % ******************* GET_OPTIMAL_PARAMETERS ******************* %
        function p = get_optimal_parameters(this, lb, ub)
        % GET_OPTIMAL_PARAMETERS returns a parameter estimation resulting from the weighted least square method applied to the accomplished measurements.
        %
        % Example:
        %     P = SOLVER_OBJECT.GET_OPTIMAL_PARAMETERS(LB, UB)
        %
        % Input:
        %     LB: the lower bound of the model parameters (optional, 
        %         default: - Inf)
        %     UB: the upper bound of the model parameters (optional,
        %         default: + Inf)
        %
        % Output:
        %     P: returns a parameter estimation resulting from the weighted
        %        least square method applied to the accomplished measurements,
        %        so that the solution is in the range LB <= P <= UB
        %
        % Throws:
        %     an error if the number of accomplished measurements is less
        %     than the dimension of the model parameters.
        %
        % The model, an initial estimation of the parameter and the
        % accomplished measurements must have been set via the associated
        % SET methods. 
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION and SET_ACCOMPLISHED_MEASUREMENTS
        %
            
            % check lower and upper bounds
            if nargin < 2 || all(lb == -Inf)
                lb = [];
            end
            
            if nargin < 3 || all(ub == Inf)
                ub = [];
            end
            
            p0 = this.get_p0();
            
            if not(isempty(lb))
                lb = util.make_column_vector(lb);
                if not(all(size(p0) == size(lb)))
                    error(this.get_message_identifier('get_optimal_parameters', 'wrong_lower_bound_size'), ['The size of the lower bound is ', mat2str(size(lb)), ' but is has to be ', mat2str(size(p0)), '.']);
                end
            end
            if not(isempty(ub))
                ub = util.make_column_vector(ub);
                if not(all(size(p0) == size(ub)))
                    error(this.get_message_identifier('get_optimal_parameters', 'wrong_upper_bound_size'), ['The size of the upper bound is ', mat2str(size(ub)), ' but is has to be ', mat2str(size(p0)), '.']);
                end
            end
            
            % return optimized parameters
            if not(isequal(this.p_lb, lb)) || not(isequal(this.p_ub, ub)) || isempty(this.p)
                this.p_lb = lb;
                this.p_ub = ub;
                
                eta = this.get_eta();
                
                if length(eta) < length(p0);
                    error(this.get_message_identifier('get_optimal_parameters', 'not_enough_accomplished_measurements'), ['The number of accomplished measurements is not enough. You need at least ', int2str(length(p0)), ' accomplished measurements, but you have only ', int2str(length(eta)), '.']);
                end
                    
                model = this.get_model();
                t = this.get_t_fix();
                v = this.get_v_fix();

                n = length(t);

                v_average = sum(v) / n;
                
                solver_po_options = this.options.get_solver_po_options();
                
                tol_x = 10^-6;                
                if solver_po_options.scale_parameters()                    
                    tol_x_min = min(abs(p0));
                    tol_x_min = 10^(floor(log10(tol_x_min)) - 4);
                    tol_x_min = max([tol_x_min; eps]);
                    if tol_x_min < tol_x
                        tol_x = tol_x_min;
                    end
                end
                
                tol_fun = 10^-8;

                opt = optimset('Display', 'off', 'Jacobian','on', 'TolX', tol_x, 'TolFun', tol_fun, 'MaxFunEvals', length(p0) * 10^3, 'MaxIter', length(p0) * 10^3);
                opt.Algorithm = solver_po_options.get_solver_algorithm();
                opt.MaxIter = solver_po_options.get_max_iter();
                opt.MaxFunEvals = solver_po_options.get_max_fun_evals();
                %{
                if lb == -Inf && ub == Inf
                    opt.Algorithm = 'levenberg-marquardt';
                    opt.ScaleProblem = 'Jacobian';
                end
                %}
                if this.use_debug()
                    opt.Display = 'iter-detailed';
                    opt.Diagnostics = 'on';
                end                
                if this.use_debug(2)
                    opt.DerivativeCheck = 'on';
                end
%                 if solver_po_options.use_algorithm_levenberg_marquardt() && solver_po_options.scale_parameters()
%                     opt.ScaleProblem = 'Jacobian';
%                 end

                
                notify(this, 'event_t_changed'); % workaround for using get_dp_Mt in fu
                if solver_po_options.scale_parameters()                    
                    S = this.get_S();
                    p0_scaled = S * p0;
                    if ~ isempty(lb)
                        lb_scaled = S * lb;
                    else
                        lb_scaled = [];
                    end
                    if ~ isempty(ub)
                        ub_scaled = S * ub;
                    else
                        ub_scaled = [];
                    end
                    
                    if solver_po_options.use_algorithm_trust_region_reflective()
                        p = S^-1 * lsqnonlin(@(p) fun(p, S), p0_scaled, lb_scaled, ub_scaled, opt);
                    else
                        p = lsqnonlin(@(p) fun(p, [], lb_scaled, ub_scaled), p0_scaled, [], [], opt);
                    end
                else
                    if solver_po_options.use_algorithm_trust_region_reflective()
                        p = lsqnonlin(@(p) fun(p), p0, lb, ub, opt);
                    else
                        p = lsqnonlin(@(p) fun(p, [], lb, ub), p0, [], [], opt);
                    end
                end
                notify(this, 'event_t_changed'); % workaround for using get_dp_Mt in fun
                
                if not(isequal(this.p, p))
                    this.p = p;
                    if this.options.use_parameter_estimation()
                        notify(this, 'event_p_changed');
                    end
                end                    
            else
                p = this.p;
            end            
            
            function [F, J] = fun(p, scale_matrix, lb, ub)
                if nargin < 2 || isempty(scale_matrix)
                    scale_matrix = eye(length(p));
                end
                inverse_scale_matrix = scale_matrix^-1;
                
                if nargin < 3
                    lb = -Inf(size(p));
                end
                if nargin < 4
                    ub = Inf(size(p));
                end
                
                F =  (v_average * v.^(-1)) .* (this.get_Mt(model, inverse_scale_matrix * p, t) - eta);
                
                violation_pattern = [];
                violation = [];
                j = 1;
                if ~ isempty(lb)
                    violation_pattern{j} = p < lb;
                    violation{j} = lb - p;
                    j = j + 1;
                end
                if ~ isempty(ub)
                    violation_pattern{j} = p > ub;
                    violation{j} = p - ub;
                end
                
                penalty_factor = 10^8;
                for j = 1:length(violation)
                    if any(violation_pattern{j})
                        F = real(F) + penalty_factor * norm(violation{j}(violation_pattern{j}))^2;
                    end
                end
                
                if nargout > 1 
                    J = diag(v_average * v.^(-1)) * this.get_dp_Mt(model, inverse_scale_matrix * p, t) * inverse_scale_matrix;
                    for i = 1:length(violation)
                        if any(violation_pattern{i})
                            J_penalty = zeros(size(J));
                            for j = 1:length(p)            
                                J_penalty(:, j) = J_penalty(:, j) + violation{i}(j);
                            end
                            J = real(J) + 2 * penalty_factor * J_penalty;
                        end
                    end
                end
            end
        end
        
        % ******************* GET_OPTIMAL_MEASUREMENTS ******************* %
        function [t_opt, w_int, w_real] = get_optimal_measurements(this, A_ineq, b_ineq)
        % GET_OPTIMAL_MEASUREMENTS returns the optimal subset of the selectable measurements.
        %
        % Example:
        %     [T_OPT, W_INT, W_REAL] = SOLVER_OBJECT.GET_OPTIMAL_MEASUREMENTS(MAX)
        %
        % Input:
        %     MAX: the maximum number of measurements allowed (optional, 
        %          default: Inf)
        %
        % Output:
        %     T_OPT: the optimal subset of the selectable measurements with
        %            less than or equal MAX measurements
        %     W_INT: the optimal integer weights of the selectable
        %            measurements which follow the constraint sum(W) <= MAX,
        %            if the direct solver is chosen or an approximation otherwise.
        %     W_REAL: the optimal real weights of the selectable 
        %             measurements which follow the constraint sum(W) <= MAX,
        %             if the direct solver is chosen or [] otherwise.
        %
        %
        % Example:
        %     [T_OPT, W_INT, W_REAL] = SOLVER_OBJECT.GET_OPTIMAL_MEASUREMENTS(A, B)
        %
        % Input:
        %     A: matrix belonging to the constraint A*W <= B
        %     B: vector belonging to the constraint A*W <= B
        %
        % Output:
        %     T_OPT: the optimal subset of the selectable measurements whose
        %             optimal real weights follow the constraint A*W <= B 
        %     W_INT: the optimal integer weights of the selectable 
        %            measurements which follow (approximately) the constraint A*W <= B,
        %            if the direct solver is chosen or an approximation otherwise.
        %     W_REAL: the optimal real weights of the selectable 
        %             measurements which follow the constraint A*W <= B,
        %             if the direct solver is chosen or [] otherwise.
        %
        % The model, an initial estimation of the parameter and the
        % selectable measurements must have been set via the associated
        % SET methods. 
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION, SET_SELECTABLE_MEASUREMENTS and SET_ACCOMPLISHED_MEASUREMENTS
        %
            switch nargin
                case 1
                    [w_int, w_real] = this.get_optimal_weights();
                case 2
                    [w_int, w_real] = this.get_optimal_weights(A_ineq);
                case 3
                    [w_int, w_real] = this.get_optimal_weights(A_ineq, b_ineq);
                otherwise
                    error(this.get_message_identifier('get_optimal_measurements', 'Too_many_input_arguments'), ['The number of input arguments exceeds 2.']);
            end
            
            t_var = this.get_t_var();
            t_opt = t_var(w_int, :);
        end                
        
        % ******************* GET_OPTIMAL_WEIGHTS ******************* %
        function [w_int, w_real] = get_optimal_weights(this, A_ineq, b_ineq)
        % GET_OPTIMAL_WEIGHTS returns the optimal weights of the selectable measurements.
        %
        % Example:
        %     [W_INT, W_REAL] = SOLVER_OBJECT.GET_OPTIMAL_WEIGHTS(MAX)
        %
        % Input:
        %     MAX: maximum number of measurements allowed (optional, 
        %         default: Inf)
        %
        % Output:
        %     W_INT: the optimal integer weights of the selectable measurements
        %            which follow the constraint sum(W) <= MAX,
        %            if the direct solver is chosen or an approximation otherwise.
        %     W_REAL: the optimal real weights of the selectable measurements
        %             which follow the constraint sum(W) <= MAX,
        %             if the direct solver is chosen or [] otherwise.
        %
        %
        % Example:
        %     [W_INT, W_REAL] = SOLVER_OBJECT.GET_OPTIMAL_WEIGHTS(A, B)
        %
        % Input:
        %     A: matrix belonging to the constraint A*W <= B
        %     B: vector belonging to the constraint A*W <= B
        %
        % Output:
        %     W_INT: the optimal integer weights of the selectable measurements
        %            which follow (approximately) the constraint A*W <= B,
        %            if the direct solver is chosen or an approximation otherwise.
        %     W_REAL: the optimal real weights of the selectable measurements
        %             which follow the constraint A*W <= B,
        %             if the direct solver is chosen or [] otherwise.
        %
        % The model, an initial estimation of the parameter and the
        % selectable measurements must have been set via the associated
        % SET methods. 
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION, SET_SELECTABLE_MEASUREMENTS and SET_ACCOMPLISHED_MEASUREMENTS
        %
            
            dim_var = length(this.t_var);
            
            % max passed
            if nargin <= 2
                if nargin == 2
                    max = A_ineq;
                else
                    max = dim_var;
                end
                if ~isscalar(max)
                    error(this.get_message_identifier('get_optimal_weights', 'max_is_no_scalar'), 'The input parameter "max" is not an scalar.');
                end
                if max > dim_var
                    warning(this.get_message_identifier('get_optimal_weights', 'wrong_max_passed'), ['max can not exceed the number of weights. max has been set to ' num2str(dim_var)']);
                    max = dim_var;
                end
                
                A_ineq = [];
                b_ineq = [];
                
                A_eq(1,:) = ones(1, dim_var);
                b_eq = max;
                
                w0 = (max / dim_var) * ones(dim_var, 1);
                
            % A, b passed    
            elseif nargin == 3
                A_eq = [];
                b_eq = [];
                w0 = ones(dim_var, 1);
            else
                error(this.get_message_identifier('get_optimal_weights', 'wrong_number_of_input_parameters'), 'You have passed the wrong number of input parameters to this method.');
            end
                
            solver_edo_options = this.options.get_solver_edo_options();
            solver_algorithm = solver_edo_options.get_solver_algorithm();
            switch solver_algorithm
                
                % local sqp solver for relaxed problem
                %case {solver_edo_options.algorithm_local_sqp, solver_edo_options.algorithm_global_sqp_serial, solver_edo_options.algortihm_global_sqp_parallel}
                case solver_edo_options.algorithm_local_sqp
                    
                    % create optimization problem structure
                    lb = zeros(dim_var, 1);
                    ub = ones(dim_var, 1);

                    tol_x = 10^-4;
                    opt = optimset('Algorithm', 'sqp', 'GradObj', 'on', 'TolCon', 0, 'TolX', tol_x);
                    opt.MaxIter = solver_edo_options.get_max_iter();
                    opt.MaxFunEvals = solver_edo_options.get_max_fun_evals();

                    if this.use_debug()
                        opt.Display = 'iter-detailed';
                        opt.Diagnostics = 'on';
                    else
                        opt.Display = 'off';
                    end
                    
                    if this.use_debug(2)
                        opt.DerivativeCheck = 'on';
                    end
                    
                    w_real = fmincon(@objective_function, w0, A_ineq, b_ineq, A_eq, b_eq, lb, ub, [], opt);

                    % calculate integer solution from real solution
                    w_int = zeros(length(w_real), 1);            
                    number_of_measurements = round(sum(w_real));
                    [~, I] = sort(w_real, 'descend');
                    for i=1:number_of_measurements
                        w_int(I(i)) = 1;
                    end
                    
                % solver direct
                case solver_edo_options.algorithm_direct                
                    w_opt = [];
                    fw_opt = [];
                    
                    check_not_ineq = isempty(A_ineq) || isempty(b_ineq);
                    check_not_eq = isempty(A_eq) || isempty(b_eq);

                    % iterate w
                    for i = 0:2^dim_var - 1
                        w = bitget(i, 1:dim_var)';

                        % check constraints
                        if (check_not_ineq || A_ineq * w <= b_ineq) && (check_not_eq || A_eq * w == b_eq)
                            try
                                fw = objective_function(w);
                                if isempty(w_opt) || fw < fw_opt
                                    w_opt = w;
                                    fw_opt = fw;
                                end
                            catch exception
                                warning(this.get_message_identifier('get_optimal_weights', 'model_error'), ['The following error occurred at w = ', num2str(w), ': ']);
                                display(getReport(exception, 'extended'));
                            end
                        end
                    end

                    w_int = w_opt;
                    w_real = [];
                    
                otherwise
                    error(this.get_message_identifier('get_optimal_weights', 'solver_unsupported'), ['The solver ' solver_edo_options.get_solver_algorithm() ' is not supported.'])
                    
            end
            
            w_int = logical(w_int);
             
            
            function [fw, dw_fw] = objective_function(w_var)
               
                if nargin >= 1
                    this.set_w_var(w_var);
                end

                % calculate function value
                model = this.get_model();
                p = this.get_p();
                t = this.get_t();
                v = this.get_v();
                w_var = this.get_w_var();
                if this.options.get_solver_edo_options().scale_parameters()
                    S = this.get_S();
                else
                    S = eye(length(p));
                end
                criterion = this.get_criterion();
                if this.use_estimation_method_region()
                    g = this.get_g();
                else
                    g = [];
                end

                fw = this.get_phiR(model, p, t, v, w_var, S, criterion, g);

                % calculate the derivative if necessary
                if nargout >= 2
                    dw_fw = this.get_dw_phiR(model, p, t, v, w_var, S, criterion, g);
                end
                
            end
            
        end
        
        
        % ******************* GET_QUALITY ******************* %
        function quality = get_quality(this, w_var)
        % GET_QUALITY returns the quality resulting from the accomplished and the by the passed weights W_VAR selected measurements.
        % The smaller the value, the better the quality.
        %
        % Example:
        %     QUALITY = SOLVER_OBJECT.GET_QUALITY(W_VAR)
        %
        % Input:
        %     W_VAR: the weights of the measurements (optional, 
        %            default: zeroes)
        %
        % Output:
        %     QUALITY: the quality resulting from the accomplished and the
        %              by the passed weights W_VAR selected measurements
        %
        % The model, an initial estimation of the parameter and the
        % selectable measurements must have been set via the associated
        % SET methods. 
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION, SET_SELECTABLE_MEASUREMENTS and SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if nargin < 2
                dim_var = length(this.t_var);
                w_var = zeros(dim_var, 1);
            end
            
            this.set_w_var(w_var);
            w_var = this.get_w_var();
            
            model = this.get_model();
            p = this.get_p();
            t = this.get_t();
            v = this.get_v();
            if this.options.get_solver_edo_options().scale_parameters()
                S = this.get_S();
            else
                S = eye(length(p));
            end
            criterion = this.get_criterion();
            if this.use_estimation_method_region()
                g = this.get_g();
            else
                g = [];
            end
            
            quality = this.get_phiR(model, p, t, v, w_var, S, criterion, g);
            
            if ~isreal(quality)
                quality = NaN;
            end
        end
        
        % ******************* GET_COVARIANCE_MATRIX ******************* %
        function C = get_covariance_matrix(this, w_var)
        % GET_COVARIANCE_MATRIX returns the covarinace matrix of the parameters resulting from the accomplished and the by the passed weights W_VAR selected measurements.
        %
        % Example:
        %     C = SOLVER_OBJECT.GET_COVARIANCE_MATRIX(W_VAR)
        %
        % Input:
        %     W_VAR: the weights of the measurements (optional, 
        %            default: zeroes)
        %
        % Output:
        %     C: the covarinace matrix of the parameters resulting from the
        %        accomplished and the by the passed weights W_VAR selected
        %        measurements
        %
        % The model, an initial estimation of the parameter and the
        % selectable measurements must have been set via the associated
        % SET methods. 
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION, SET_SELECTABLE_MEASUREMENTS and SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if nargin < 2
                dim_var = length(this.t_var);
                w_var = zeros(dim_var, 1);
            end
            
            this.set_w_var(w_var);
            w_var = this.get_w_var();
            
            model = this.get_model();
            p = this.get_p();
            t = this.get_t();
            v = this.get_v();
            S = eye(length(p));
            
            C = this.get_C(model, p, t, v, w_var, S);
        end
        
        % ******************* GET_CONFIDENCE_INTERVALS ******************* %
        function v = get_confidence_intervals(this, w_var, alpha)
        % GET_CONFIDENCE_INTERVALS returns the absolute values by which the optimal parameters maximally vary from the true parameters with a probability of alpha.
        %
        % Example:
        %     V = SOLVER_OBJECT.GET_CONFIDENCE_INTERVALS(W_VAR, ALPHA)
        %
        % Input:
        %     W_VAR: the weights of the measurements (optional, 
        %            default: zeroes)
        %     ALPHA: the probability that the returned values ​​should be
        %            accurate (optional, default: the value which is set in 
        %            the options of the solver object)
        %
        % Output:
        %     V: the absolute values by which the optimal parameters
        %        maximally vary from the true parameters with a probability
        %        of alpha.
        %
        % The model, an initial estimation of the parameter and the
        % selectable measurements must have been set via the associated
        % SET methods. 
        %
        % see also SET_MODEL, SET_INITIAL_PARAMETER_ESTIMATION, SET_SELECTABLE_MEASUREMENTS and SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if nargin < 2 || isempty(w_var)
                C = this.get_covariance_matrix();
            else
                C = this.get_covariance_matrix(w_var);
            end
            
            if nargin < 3
                alpha = this.options.get_alpha();
            end
            g = chi2inv(alpha, length(this.get_p()));
            
            v = sqrt(diag(C) * g);
        end
        
        
        % ******************* THE SETTER FOR THE INPUTS ******************* %
        function set_model(this, model)
        % SET_MODEL sets the model which will be used for the computations.
        %
        % Example:
        %     SOLVER_OBJECT.SET_MODEL(MODEL)
        %
        % Input:
        %     MODEL: the model represented as an object whose class
        %            implements the MODEL interface
        %
        % see also MODEL, GET_OPTIMAL_WEIGHTS, GET_QUALITY and GET_OPTIMAL_PARAMETERS
        %
        
            if not(isequal(this.model, model))
                this.model = model;    
                notify(this, 'event_model_changed');
            end
        end
        
        function set_initial_parameter_estimation(this, p0)
        % SET_INITIAL_PARAMETER_ESTIMATION sets the initial estimation of the model parameters.
        %
        % Example:
        %     SOLVER_OBJECT.SET_INITIAL_PARAMETER_ESTIMATION(P0)
        %
        % Input:
        %     P0: the initial estimation of the model parameters
        %
        % see also GET_OPTIMAL_WEIGHTS, GET_QUALITY and GET_OPTIMAL_PARAMETERS
        %
        
            p0 = util.make_column_vector(p0);            
            if not(isequal(this.p0, p0))
                this.p0 = p0;    
                notify(this, 'event_p0_changed');
            end
        end
        
        function set_selectable_measurements(this, t_var, v_var)
        % SET_SELECTABLE_MEASUREMENTS sets the selectable measurements and the associated variances of the measurement errors.
        % These are the measurements for which the quality will be
        % calculated or optimized. 
        %
        % Example:
        %     SOLVER_OBJECT.SET_SELECTABLE_MEASUREMENTS(T_VAR, V_VAR)
        %
        % Input:
        %     T_VAR: the experimental designs for which the quality will be
        %            calculated or optimized
        %            format: a n x m matrix where n is the number of
        %                    measurements and m the dimension of a
        %                    measurement
        %     V_VAR: the variances of the associated measurement errors
        %            format: a vector of length n where n is the number of 
        %                    measurements
        %
        % see also GET_OPTIMAL_WEIGHTS, GET_QUALITY and GET_OPTIMAL_PARAMETERS
        %    
        
            v_var = util.make_column_vector(v_var);
            size_t_var = size(t_var); 
            if size_t_var(1) ~= length(v_var)
                error(this.get_message_identifier('set_selectable_measurements', 'dimension_mismatch'), 'The dimensions of t_var and v_var do not match.');
            end                
            
            if not(isequal(this.t_var, t_var))
                this.t_var = t_var;
                notify(this, 'event_t_changed');
            end
            
            if not(isequal(this.v_var, v_var))
                this.v_var = v_var;               
                notify(this, 'event_v_changed');
            end
        end
        
        function set_accomplished_measurements(this, t_fix, v_fix, eta)
        % SET_ACCOMPLISHED_MEASUREMENTS sets the accomplished measurements including the measurements and the associated variances of the measurement errors.
        %
        % Example:
        %     SOLVER_OBJECT.SET_ACCOMPLISHED_MEASUREMENTS(T_FIX, V_FIX, ETA)
        %
        % Input:
        %     T_FIX: the experimental designs of accomplished measurements
        %            format: a n x m matrix where n is the number of
        %                    measurements and m the dimension of a
        %                    measurement
        %     V_FIX: the variances of the associated measurement errors
        %            format: a vector of length n where n is the number of 
        %                    measurements
        %     ETA: the associated measurement results (optional,
        %          only necessary for parameter estimation)
        %          format: a vector of length n where n is the number of 
        %                  measurements
        %
        % see also GET_OPTIMAL_WEIGHTS, GET_QUALITY and GET_OPTIMAL_PARAMETERS
        %
        
            v_fix = util.make_column_vector(v_fix);
            size_t_fix = size(t_fix); 
            
            if size_t_fix(1) ~= length(v_fix)
                error(this.get_message_identifier('set_accomplished_measurements', 'dimension_mismatch'), 'The dimensions of t_fix and v_fix do not match.');
            end                
            
            if nargin >= 4
                eta = util.make_column_vector(eta);

                if length(v_fix) ~= length(eta)
                    error(this.get_message_identifier('set_accomplished_measurements', 'dimension_mismatch'), 'v_fix and eta must have the same dimensions.');
                end
            else
                eta = [];
            end
            
            
            if not(isequal(this.t_fix, t_fix))
                this.t_fix = t_fix;
                notify(this, 'event_t_changed');
            end
            
            if not(isequal(this.v_fix, v_fix))
                this.v_fix = v_fix;
                notify(this, 'event_v_changed');
            end

            if not(isequal(this.eta, eta))
                this.eta = eta;
                notify(this, 'event_eta_changed');
            end
        end
        
        
        function set_option(this, name, value)
        % SET_OPTION sets the option for this SOLVER object.
        %
        % Example:
        %     SOLVER_OBJECT.SET_OPTION(NAME, VALUE)
        %
        % Input:
        %     NAME: the name of the option to be changed
        %     VALUE: the new value of the option
        %
        %     Possible options are:
        %     'parameter_estimation': whether a parameter estimation should
        %         be performed before the optimal design estimation or not
        %         (possible values: 'yes', 'no', default: 'no')
        %     'estimation_method': the method of the estimation of the
        %         quality of the experimental design (possible values: 
        %         'point', 'region', default: 'region')
        %     'alpha': the confidence level for the region estimation 
        %         (possible value: a scalar with 0 < alpha < 1,
        %         default: 0.95)
        %     'debug': the level of debug information to be displayed
        %         (possible value: a non-negative integer,
        %         default: 0 (no debug informations))
        %     'criterion': the criterion for the quality of an experimental
        %         design (possible values: an object of the CRITERION class,
        %         default: a CRITERION_A object)
        %     'solver_edo_options': the solver to be used to solve the
        %         experimental design optimization problem (possible
        %         values: an object of the SOLVER_EDO_OPTIONS class,
        %         default: a default SOLVER_EDO_OPTIONS object)
        %     'edo_algorithm': the method to be used to solve the experimental
        %         design optimization problem (possible values: 'direct', 
        %         'local_sqp', default: 'local_sqp')
        %     'edo_scale_parameters': whether the parameters should be scaled 
        %         or not (possible values: 'yes', 'no', default: 'yes')
        %     'edo_max_fun_evals': the number of maximal model evaluations done
        %         by the 'local_sqp' solver (possible values: a non-negative integer,
        %         default: 10^3)
        %     'edo_max_iter': the number of maximal iterations of the 'local_sqp' solver
        %         (possible values: a  non-negative integer, default: 10^3)
        %     'po_algorithm': the method to be used to solve the parameter
        %         optimization problem (possible values: 'trust-region-reflectiv', 
        %         'levenberg-marquardt', default: 'trust-region-reflectiv')
        %     'po_scale_parameters': whether the parameters have to be scaled
        %         for the optimization (possible values: 'yes', 'no',
        %         default: 'no')
        %     'po_max_fun_evals': the number of maximal model evaluations done
        %         by the solver (possible values: a non-negative integer,
        %         default: 3 * 10^3)
        %     'po_max_iter': the number of maximal iterations of the solver
        %         (possible values: a  non-negative integer, default: 3 * 10^3)
        %
        % Throws:
        %     An error if a value doesn't match to an option or a wrong
        %     option is passed.
        %
           
            event_id = [];
            switch name
                case solver_options.alpha_id
                    if ~ isequal(this.options.get_alpha(), value)
                        event_id = 'event_alpha_changed';
                    end
                case solver_options.criterion_id
                    if ~ isequal(this.options.get_criterion(), value)
                        event_id = 'event_criterion_changed';
                    end
                case solver_options.estimation_method_id
                    if ~ isequal(this.options.get_estimation_method(), value)
                        event_id = 'event_estimation_method_changed';
                    end
                case solver_options.parameter_estimation_id
                    if ~ isequal(this.options.use_parameter_estimation(), value)
                        event_id = 'event_p_changed';
                    end
                case solver_options.get_solver_edo_options().scale_parameters_id
                    if ~ isequal(this.options.get_solver_edo_options().scale_parameters(), value)
                        event_id = 'event_v_changed';
                    end
            end
            if strncmp(name, 'po_', 3)
                event_id = 'event_p_changed';
            end
           
            this.options.set_option(name, value);
            if ~ isempty(event_id)
                notify(this, event_id)
            end
        end
                
        
        function g = get_g(this)
        % GET_G returns the radius of the neighborhood of the parameters P for which the weights of the measurements will be optimized.
        %
        % Example:
        %     G = SOLVER_OBJECT.GET_G()
        %
        % Output:
        %     G: the radius of the neighborhood of the parameters P for
        %        which the weights of the measurements will be
        %        optimized
        %
        % Throws:
        %     an error if g is not set.
        %
        % see also SET_ALPHA
        %
        
            alpha = this.options.get_alpha();
            g = chi2inv(alpha, length(this.get_p()));
        end
         
    end
    
    
    methods (Access = protected)
                
        
        % ******************* INTERNAL CALCULATION METHODS ******************* %
        
        function Mt = get_Mt(this, model, p, t)
        % GET_MT returns the output of the MODEL with the parameters P at the measurements T.
        %
        % Example:
        %     MT = SOLVER_OBJECT.GET_MT(MODEL, P, T)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements for which the output will be 
        %        calculated
        %
        % Output:
        %     MT: the output of the MODEL with the parameters P at the
        %         measurements T.
        %
        
            size_t = size(t);
            m = size_t(1);
            Mt = zeros(m, 1);
            for i=1:m
                Mt(i) = model.get_M(p, this.get_ti(t, i));
            end  
            
        end
                
        function dp_Mt = get_dp_Mt(this, model, p, t)
        % GET_DP_MT returns the first derivation of the MODEL with respect to the parameters P at the measurements T.
        %
        % Example:
        %     DP_MT = SOLVER_OBJECT.GET_DP_MT(MODEL, P, T)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements for which the derivations will be 
        %        calculated
        %
        % Output:
        %     DP_MT: the first derivation of the MODEL with respect to the 
        %            parameters P at the measurements T
        %
        
            if isempty(this.dp_Mt)
                dim_p = length(p);
                size_t = size(t);
                m = size_t(1);
                dp_Mt = zeros(m, dim_p);

                for i=1:m
                    dp_Mt(i,:) = model.get_dp_M(p, this.get_ti(t, i))';
                end        

                this.dp_Mt = dp_Mt;
            else
                dp_Mt = this.dp_Mt;
            end
        end                        
                
        function dpdp_Mt = get_dpdp_Mt(this, model, p, t)
        % GET_DPDP_MT returns the second derivation of the MODEL with respect to the parameters P at the measurements T.
        %
        % Example:
        %     DPDP_MT = SOLVER_OBJECT.GET_DPDP_MT(MODEL, P, T)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements for which the derivations will be 
        %        calculated
        %
        % Output:
        %     DPDP_MT: the second derivation of the MODEL with respect to 
        %              the parameters P at the measurements T
        %
        
            if isempty(this.dpdp_Mt)
                dim_p = length(p);
                size_t = size(t);
                m = size_t(1);
                dpdp_Mt = cell(dim_p, 1);

                for i=1:dim_p
                    dpdp_Mt{i} = zeros(m, dim_p);
                end

                for i=1:m
                    A = model.get_dpdp_M(p, this.get_ti(t, i));

                    if not(isempty(A)) && all(all(isfinite(A)))
                        for j=1:dim_p
                            dpdp_Mt{j}(i,:) = A(j,:);
                        end
                    else
                        if isempty(A)
                            error(this.get_message_identifier('get_dpdp_Mt', 'empty_matrix'), 'model.get_dpdp_M returned an empty matrix.')
                        else
                            error(this.get_message_identifier('get_dpdp_Mt', 'NaN'), 'model.get_dpdp_M returned NaN.')
                        end
                    end                                         
                end

                this.dpdp_Mt = dpdp_Mt;
            else
                dpdp_Mt = this.dpdp_Mt;
            end
        end        
        
        
        function C = get_C(this, model, p, t, v, w_var, S)
        % GET_C returns the associated covariance matrix.
        %
        % Example:
        %     C = SOLVER_OBJECT.GET_C(MODEL, P, T, V, W)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %
        % Output:
        %     C: the associated covariance matrix
        %
        
            if isempty(this.C)
                dp_Mt = this.get_dp_Mt(model, p, t);
                
                dim_fix = length(t) - length(w_var);
                w = [ones(dim_fix, 1); w_var];                
                dim_w = length(w);
                VW = spdiags(w .* v.^(-1), 0, dim_w, dim_w);
                C_inv = dp_Mt' * VW * dp_Mt;
                warning('off', 'MATLAB:singularMatrix');
                C = inv(C_inv);
                warning('on', 'MATLAB:singularMatrix');
                if not(any(any(isinf(C))))
                    C = S * C * S;
                else
                    C = ones(size(C)) * inf;
                end
                this.C = C;
            else
                C = this.C;
            end
        end
                
        function dw_C = get_dw_C(this, model, p, t, v, w_var, S)
        % GET_DW_C returns the derivation of the associated covariance matrix with respect to the weights W.
        %
        % Example:
        %     DW_C = SOLVER_OBJECT.GET_DW_C(MODEL, P, T, V, W)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %
        % Output:
        %     DW_C: the derivation of the associated covariance matrix with
        %           respect to the weights W
        %
        
            if isempty(this.dw_C)
                dp_Mt = this.get_dp_Mt(model, p, t); 
                C = this.get_C(model, p, t, v, w_var, S);
                A = dp_Mt * S^-1 * C;

                dim_w_var = length(w_var);
                dim_fix = length(t) - dim_w_var;            
                dw_C = cell(dim_w_var, 1);

                for i = 1:dim_w_var
                    j = i + dim_fix;
                    dw_C{i} = - A(j,:)' * v(j)^-1 * A(j,:) ;
                end

                this.dw_C = dw_C;
            else
                dw_C = this.dw_C;
            end
        end
                
        function dp_C = get_dp_C(this, model, p, t, v, w_var, S)
        % GET_DP_C returns the derivation of the associated covariance matrix with respect to the parameters P.
        %
        % Example:
        %     DP_C = SOLVER_OBJECT.GET_DP_C(MODEL, P, T, V, W)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %
        % Output:
        %     C: the derivation of the associated covariance matrix with 
        %        respect to the parameters P
        %
        
            if isempty(this.dp_C)
                C = this.get_C(model, p, t, v, w_var, S);
                H = this.get_H(model, p, t, v, w_var, S);
                
                dim_p = length(p);
                dp_C = cell(dim_p, 1);
                
                for i=1:dim_p
                    dp_C{i} = - C * (H{i}' + H{i}) * C;
                end
                
                this.dp_C = dp_C;
            else
                dp_C = this.dp_C;
            end
        end
                        
        
        function H = get_H(this, model, p, t, v, w_var, S)
        % GET_H returns the associated auxiliary matrix H.
        %
        % Example:
        %     H = SOLVER_OBJECT.GET_H(MODEL, P, T, V, W)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %
        % Output:
        %     H: the auxiliary matrix H
        %
        
            if isempty(this.H)
                dp_Mt = this.get_dp_Mt(model, p, t);
                dpdp_Mt = this.get_dpdp_Mt(model, p, t);
                
                dim_p = length(p);
                H = cell(dim_p, 1);
                
                dim_fix = length(t) - length(w_var);
                w = [ones(dim_fix, 1); w_var];
                dim_w = length(w);
                VW = spdiags(w .* v.^(-1), 0, dim_w, dim_w);
                
                G = S^-1 * dp_Mt' * VW;
                
                for i=1:dim_p
                    H{i} =  G * dpdp_Mt{i} * S^-1;
                end
                this.H = H;
            else
                H = this.H;
            end
        end        
        
        function dw_H = get_dw_H(this, model, p, t, v, w_var, S)
        % GET_H returns the derivation of the associated auxiliary matrix H with respect to the weights W.
        %
        % Example:
        %     DW_H = SOLVER_OBJECT.GET_DW_H(MODEL, P, T, V, W)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %
        % Output:
        %     dw_H: the derivation of the associated auxiliary matrix H 
        %           with respect to the weights W
        %
            
            if isempty(this.dw_H)
                dp_Mt = this.get_dp_Mt(model, p, t);
                dpdp_Mt = this.get_dpdp_Mt(model, p, t);

                dim_p = length(p);
                dim_w_var = length(w_var);
                dim_fix = length(t) - dim_w_var; 
                dw_H = cell(dim_p, dim_w_var);                 

                A = S^-1 * dp_Mt';
                for i = 1:dim_p
                    B = dpdp_Mt{i} * S^-1;
                    for j = 1:dim_w_var
                        k = j + dim_fix;
                        dw_H{i,j} =  A(:,k) * v(k)^(-1) * B(k,:); 
                    end
                end

                this.dw_H = dw_H;
            else
                dw_H = this.dw_H;
            end
        end
        
        
        function phi = get_phi(this, model, p, t, v, w_var, S, criterion)
        % GET_PHI returns the quality resulting of the weights W, for a model in which the parameters appear linearly.
        %
        % Example:
        %     PHI = SOLVER_OBJECT.GET_PHI(MODEL, P, T, V, W, CRITERION)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %     CRITERION: a class that implements the CRITERION interface
        %
        % Output:
        %     PHI: the quality resulting of the weights W
        %
        
            if isempty(this.phi)
                C = this.get_C(model, p, t, v, w_var, S);
                phi = criterion.get_phi(C);
                this.phi = phi;
            else
                phi = this.phi;
            end
        end    
        
        function dp_phi = get_dp_phi(this, model, p, t, v, w_var, S, criterion)
        % GET_DP_PHI returns the derivation with respect to the parameters P of the quality resulting of the weights W, for a model in which the parameters appear linearly.
        %
        % Example:
        %     DP_PHI = SOLVER_OBJECT.GET_DP_PHI(MODEL, P, T, V, W, CRITERION)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %     CRITERION: a class that implements the CRITERION interface
        %
        % Output:
        %     DP_PHI: the derivation with respect to the parameters P of
        %             the quality resulting of the weights W
        %
        
            if isempty(this.dp_phi)
                C = this.get_C(model, p, t, v, w_var, S);
                dp_C = this.get_dp_C(model, p, t, v, w_var, S);

                dim_p = length(p);
                dp_phi = zeros(dim_p, 1);

                for i=1:dim_p
                    dp_phi(i) = criterion.get_dCd_phi(C, dp_C{i});
                end
                this.dp_phi = dp_phi;
            else
                dp_phi = this.dp_phi;
            end
        end        
                
        function dw_phi = get_dw_phi(this, model, p, t, v, w_var, S, criterion)
        % GET_DW_PHI returns the derivation with respect to the weights W of the quality resulting of the weights W, for a model in which the parameters appear linearly.
        %
        % Example:
        %     DW_PHI = SOLVER_OBJECT.GET_DW_PHI(MODEL, P, T, V, W, CRITERION)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %     CRITERION: a class that implements the CRITERION interface
        %
        % Output:
        %     DW_PHI: the derivation with respect to the weights W of 
        %             the quality resulting of the weights W
        %
        
            if isempty(this.dw_phi)
                C = this.get_C(model, p, t, v, w_var, S);
                dw_C = this.get_dw_C(model, p, t, v, w_var, S);

                dim_w_var = length(w_var);
                dw_phi = zeros(dim_w_var, 1);

                for i = 1:dim_w_var
                    dw_phi(i) = criterion.get_dCd_phi(C, dw_C{i});
                end
                this.dw_phi = dw_phi;
            else
                dw_phi = this.dw_phi;
            end
        end
                
        function dwdp_phi = get_dwdp_phi(this, model, p, t, v, w_var, S, criterion)
        % GET_DWDP_PHI returns the derivation with respect to the parameters P and then to the weights W of the quality resulting of the weights W, for a model in which the parameters appear linearly.
        %
        % Example:
        %     DWDP_PHI = SOLVER_OBJECT.GET_DWDP_PHI(MODEL, P, T, V, W, CRITERION)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %     CRITERION: a class that implements the CRITERION interface
        %
        % Output:
        %     DWDP_PHI: the derivation with respect to the parameters P and
        %               then to the weights W of the quality resulting
        %               of the weights W
        %
        
            if isempty(this.dwdp_phi)
                C = this.get_C(model, p, t, v, w_var, S);
                dp_C = this.get_dp_C(model, p, t, v, w_var, S);
                dw_C = this.get_dw_C(model, p, t, v, w_var, S);
                H = this.get_H(model, p, t, v, w_var, S);
                dw_H = this.get_dw_H(model, p, t, v, w_var, S);

                dim_p = length(p);
                dim_w_var = length(w_var);
                dwdp_phi = zeros(dim_p, dim_w_var);

                for i = 1:dim_p
                    for j = 1:dim_w_var
                        dw_dp_C_i_j = - (dw_C{j} * (H{i}' + H{i}) * C + C * (dw_H{i,j}' + dw_H{i,j}) * C + C * (H{i}' + H{i}) * dw_C{j});
                        dwdp_phi(i,j) = criterion.get_dCd_dCd_phi(C, dp_C{i}, dw_C{j}) + criterion.get_dDd_dCd_phi(C, dp_C{i}, dw_dp_C_i_j);
                    end
                end
                this.dwdp_phi = dwdp_phi;
            else
                dwdp_phi = this.dwdp_phi;
            end
        end
        
        
        function C0 = get_C0(this, model, p, t, v, w_var, S)
        % GET_C0 returns the associated covariance matrix with w_var = 0.
        %
        % Example:
        %     C0 = SOLVER_OBJECT.GET_C0(MODEL, P, T, V, W_VAR, S)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %
        % Output:
        %     C0: the associated covariance matrix with w_var = 0
        %
        
            if isempty(this.C0)
                C0 = this.get_C(model, p, t, v, zeros(length(w_var), 1), S);
                this.C0 = C0;
            else
                C0 = this.C0;
            end
        end
        
        
        function phiR = get_phiR(this, model, p, t, v, w_var, S, criterion, g)
        % GET_PHIR returns the quality resulting of the weights W.
        %
        % Example:
        %     PHIR = SOLVER_OBJECT.GET_PHIR(MODEL, P, T, V, W, S, CRITERION, G)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %     CRITERION: a class that implements the CRITERION interface
        %     g: the radius of the neighborhood of the parameters P for
        %        which the weights of the measurements will be
        %        optimized
        %
        % Output:
        %     PHIR: the quality resulting of the weights W
        %
        
            if isempty(this.phiR)
                phi = this.get_phi(model, p, t, v, w_var, S, criterion);
                
                if this.use_estimation_method_region()
                    C0 = this.get_C0(model, p, t, v, w_var, S);
                    dp_phi = this.get_dp_phi(model, p, t, v, w_var, S, criterion);

                    phiR = phi + g^(1/2) * (dp_phi' * C0 * dp_phi)^(1/2);
                else
                    phiR = phi;
                end
                
                this.phiR = phiR;
            else
                phiR = this.phiR;
            end
        end
        
        function dw_phiR = get_dw_phiR(this, model, p, t, v, w_var, S, criterion, g)
        % GET_DW_PHIR returns the derivation with respect to the weights W of the quality resulting of the weights W.
        %
        % Example:
        %     DW_PHIR = SOLVER_OBJECT.GET_DW_PHIR(MODEL, P, T, V, W, S, CRITERION, G)
        %
        % Input:
        %     MODEL: a class that implements the MODEL interface
        %     P: the parameters for the MODEL
        %     T: the measurements
        %     V: the variances associated with these measurements
        %     W_VAR: the weights for the variable measurements
        %     S: the scaling matrix for the parameters
        %     CRITERION: a class that implements the CRITERION interface
        %     g: the radius of the neighborhood of the parameters P for
        %        which the weights of the measurements will be
        %        optimized
        %
        % Output:
        %     DW_PHIR: the derivation with respect to the weights W of
        %              the quality resulting of the weights W
        %

            if isempty(this.dw_phiR)                     
                dw_phi = this.get_dw_phi(model, p, t, v, w_var, S, criterion);
                
                if this.use_estimation_method_region()
                    C0 = get_C0(this, model, p, t, v, w_var, S);
                    C = this.get_C(model, p, t, v, w_var, S);
                    dw_C = this.get_dw_C(model, p, t, v, w_var, S);
                    dp_phi = this.get_dp_phi(model, p, t, v, w_var, S, criterion);
                    dwdp_phi = this.get_dwdp_phi(model, p, t, v, w_var, S, criterion);                

                    m = length(dw_C);
                    dw_phiR = zeros(m, 1);

                    s = (dp_phi' * C0 * dp_phi);
                    if s > 0
                        for i=1:m
                            dw_phiR(i) = dw_phi(i) + g^(1/2) * s^(-1/2) * dwdp_phi(:, i)' * C * dp_phi;
                        end
                    else
                        dw_phiR = dw_phi;
                    end
                else
                    dw_phiR = dw_phi;
                end
                                
                this.dw_phiR = dw_phiR;
            else
                dw_phiR = this.dw_phiR;
            end
            
        end
        
        
        
        % ******************* INTERNAL GETTER AND SETTER ******************* %
        
        function model = get_model(this)
        % GET_MODEL returns the model which should be used for the computations.
        %
        % Example:
        %     MODEL = SOLVER_OBJECT.GET_MODEL()
        %
        % Output:
        %     MODEL: a class that implements the MODEL interface
        %
        % Throws:
        %     an error if the model is not set.
        %
        % see also MODEL and SET_MODEL
        %
        
            if isempty(this.model)
                error(this.get_message_identifier('get_model', 'not_set'), 'model is not set.');
            else
                model = this.model;
            end
        end
        
        function p0 = get_p0(this)
        % GET_P0 returns the initial estimation of the model parameters.
        %
        % Example:
        %     P0 = SOLVER_OBJECT.GET_P0()
        %
        % Output:
        %     P0: the initial estimation of the model parameters
        %
        % Throws:
        %     an error if the initial estimation of the model parameters
        %     are not set.
        %
        % see also SET_INITITAL_PARAMETERS
        %
        
            if isempty(this.p0)
                error(this.get_message_identifier('get_p0', 'not_set'), 'The inital estimation of the model parameters is not set.');
            else
                p0 = this.p0;
            end
        end
        
        function p = get_p(this)
        % GET_P returns the optimized parameters estimation if enough measurements are accomplished and the usage of the optimized parameter estimation is enabled in the options or the initial parameter estimation otherwise.
        %
        % Example:
        %     P = SOLVER_OBJECT.GET_P()
        %
        % Output:
        %     P: the parameters for the model
        %
        % Throws:
        %     an error if the initial estimation of the model parameters
        %     are not set.
        %
        % see also SET_INITITAL_PARAMETERS
        %
            if this.options.use_parameter_estimation()
                try
                    p = this.get_optimal_parameters();
                catch e
                    if strcmp(e.identifier, this.get_message_identifier('get_optimal_parameters', 'not_enough_accomplished_measurements'))
                        p = this.get_p0();
                    else
                        rethrow(e);
                    end
                end
            else
                p = this.get_p0();
            end
        end    
        
        function t = get_t(this)
        % GET_T returns the measurements which should be used for the computations.
        %
        % Example:
        %     T = SOLVER_OBJECT.GET_T()
        %
        % Output:
        %     T: the measurements
        %
        % Throws:
        %     an error if the variable measurements are not set.
        %
        % see also SET_T_V_VAR
        %
        
            if isempty(this.t_var)
                error(this.get_message_identifier('get_t', 'not_set'), 't_var is not set.');
            else
                t = [this.t_fix; this.t_var];
            end
        end
        
        function ti = get_ti(this, t, i)
        % GET_TI returns the i-th measurement of the passed measurements.
        %
        % Example:
        %     TI = SOLVER_OBJECT.GET_TI(T)
        %
        % Input: 
        %     T: the measurements
        %     I: the number of the desired measurement
        %
        % Output:
        %     TI = the i-th measurement
        %
        
            ti = t(i, :);
        end      
        
        function t_fix = get_t_fix(this)
        % GET_T_FIX returns the measurements from the accomplished measurements.
        %
        % Example:
        %     T_FIX = SOLVER_OBJECT.GET_T_FIX()
        %
        % Output:
        %     T_FIX: the measurements from the accomplished measurements
        %
        % Throws:
        %     an error if the measurements from the accomplished
        %     measurements are not set.
        %
        % see also SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if isempty(this.t_fix)
                error(this.get_message_identifier('get_t_fix', 'not_set'), 't_fix is not set.');
            else
                t_fix = this.t_fix;
            end
        end
        
        function t_var = get_t_var(this)
        % GET_T_VAR returns the measurements from the selectable measurements.
        %
        % Example:
        %     T_VAR = SOLVER_OBJECT.GET_T_VAR()
        %
        % Output:
        %     T_VAR: the measurements from the selectable measurements
        %
        % Throws:
        %     an error if the measurements from the selectable
        %     measurements are not set.
        %
        % see also SET_SELECTABLE_MEASUREMENTS
        %
        
            if isempty(this.t_var)
                error(this.get_message_identifier('get_t_var', 'not_set'), 't_var is not set.');
            else
                t_var = this.t_var;
            end
        end
        
        function v = get_v(this)
        % GET_V returns the variances from the measurement results.
        %
        % Example:
        %     V = SOLVER_OBJECT.GET_V()
        %
        % Output:
        %     V: the variances from the measurement results
        %
        % Throws:
        %     an error if the variances are not set.
        %
        % see also SET_T_V_VAR
        %
        
            if isempty(this.t_var)
                error(this.get_message_identifier('get_v', 'not_set'), 'v is not set.');
            else
                v = [this.v_fix; this.v_var];
            end
        end
        
        function v_fix = get_v_fix(this)
        % GET_V_FIX returns the variance from the accomplished measurements.
        %
        % Example:
        %     V_FIX = SOLVER_OBJECT.GET_V_FIX()
        %
        % Output:
        %     V_FIX: the variance from the accomplished measurements
        %
        % Throws:
        %     an error if the variance from the accomplished measurements
        %     are not set.
        %
        % see also SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if isempty(this.v_fix)
                error(this.get_message_identifier('get_v_fix', 'not_set'), 'v_fix is not set.');
            else
                v_fix = this.v_fix;
            end
        end
        
        function eta = get_eta(this)
        % GET_ETA returns the measurement results from the accomplished measurements.
        %
        % Example:
        %     ETA = SOLVER_OBJECT.GET_ETA()
        %
        % Output:
        %     ETA: the measurement results from the accomplished
        %          measurements
        %
        % Throws:
        %     an error if the measurement results from the accomplished
        %     measurements are not set.
        %
        % see also SET_ACCOMPLISHED_MEASUREMENTS
        %
        
            if isempty(this.eta)
                error(this.get_message_identifier('get_eta', 'not_set'), 'eta is not set.');
            else
                eta = this.eta;
            end
        end
        
        function w_var = get_w_var(this)
        % GET_W_VAR returns the weights of the variable measurements which should be used for the computations.
        %
        % Example:
        %     W_VAR = SOLVER_OBJECT.GET_W_VAR()
        %
        % Output:
        %     W_VAR: the weights of the variable measurements
        %
        % Throws:
        %     an error if the variable weights are not set.
        %
        % see also SET_W_VAR
        %
        
            if isempty(this.w_var)
                error(this.get_message_identifier('get_w_var', 'not_set'), 'w_var is not set.');
            else
                w_var = this.w_var;
            end
        end
        
        function set_w_var(this, w_var)
        % SET_W_VAR sets the weights for the variable measurements.
        %
        % Example:
        %     SOLVER_OBJECT.SET_W_VAR(W_VAR)
        %
        % Input:
        %     W_VAR: the weights for the variable measurements
        %
        
            w_var = util.make_column_vector(w_var);
            
            if not(all(size(w_var) == size(this.v_var)))
                error(this.get_message_identifier('set_w_var', 'not_set'), ['The weights w_var must have size ', mat2str(size(this.v_var)), ' but the size is ', mat2str(size(w_var)), '.']);                
            end
            
            w_var(w_var < 0 & w_var >= -eps) = 0;
            
            if any(w_var < 0)
                error(this.get_message_identifier('set_w_var', 'not_set'), 'One of the weights for the variable measurements is negative.');
            end
            
            if not(isequal(this.w_var, w_var))
                this.w_var = w_var;
                notify(this, 'event_w_changed');                
            end
        end
        
        function S = get_S(this)
        % GET_S returns the scaling matrix, which is the inverse of the diagonal matrix of the parameters.
        %
        % Example:
        %     S = SOLVER_OBJECT.GET_S()
        %
        % Output:
        %     S: the scaling matrix, which is the inverse of the diagonal
        %        matrix of the parameters
        %
        % Throws:
        %     an error if p is not set.
        %
        % see also SET_P
        %
        
            p = abs(this.get_p());
            S = diag(p.^(-1));
            %dim_p = length(p);
            
            %if this.options.get_solver_edo_options().scale_parameters()
            %    S = spdiags(p.^(-1), 0, dim_p, dim_p);
            %else
            %    S = spdiags(ones(dim_p,1), 0, dim_p, dim_p);
            %end
            
        end
        
        
        function boolean = use_debug(this, debug)
        % USE_DEBUG returns whether debug informations will be output or not.
        %
        % Example:
        %     SOLVER_OBJECT.USE_DEBUG(DEBUG)
        %
        % Input:
        %     DEBUG: the debug-level which will be checked
        %            (optional, default: 1)
        %
        % Output:
        %     USE_DEBUG: whether the current debug-level is greater or
        %                equal the passed debug-level
        %
        % see also SET_DEBUG
        %
        
            if nargin == 2
                boolean = this.options.use_debug(debug);
            else
                boolean = this.options.use_debug();
            end
        end
                
        
        
        function criterion = get_criterion(this)
        % GET_CRITERION returns the criterion which characterize the quality of the weights of the measurements.
        %
        % Example:
        %     CRITERION = SOLVER_OBJECT.GET_CRITERION()
        %
        % Output:
        %     CRITERION: a class that implements the CRITERION interface
        %
        % Throws:
        %     an error if the criterion is not set.
        %
        % see also CRITERION and SET_CRITERION
        %
        
            criterion = this.options.get_criterion();
        end
        
        function boolean = use_estimation_method_region(this)
        % ARE_MEASUREMENTS_ACCOMPLISHED returns whether measurements are accomplished or not.
        %
        % Example:
        %     BOOLEAN = SOLVER_OBJECT.ARE_MEASUREMENTS_ACCOMPLISHED()
        %
        % Output:
        %     BOOLEAN: whether measurements are accomplished or not
        %
        
            % sufficient measurements available && region estimation method desired
            boolean = length(this.t_fix) >= length(this.p0) && this.options.use_estimation_method_region();
        end
               
        
    end
    
    methods (Access = protected, Static)
        
        function s = get_message_identifier(method, mnemonic)
        % GET_MESSAGE_IDENTIFIER returns the identifier for an error or a warning raised in methods of these object.
        %
        % Example:
        %     ID = SOLVER.GET_MESSAGE_IDENTIFIER(METHOD, MNEMONIC)
        %
        % Input:
        %     METHOD: the method in which an error or a warning occurred
        %     MNEMONIC: a unique keyword for the error or warning
        %
        % Output:
        %     ID: the identifier for the error or a warning
        %
            s = util.get_message_identifier('solver', method, mnemonic);
        end       
        
    end
    
end

%#ok<*PROP>
