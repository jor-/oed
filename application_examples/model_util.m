classdef model_util
% MODEL_UTIL  represents an utility class with static auxiliary methods for working with model classes.
%
% MODEL_UTIL Methods:
%   GET_FICTITIOUS_MEASUREMENT_RESULTS - returns fictitious measurement results according to the passed informations.
%   GET_MODEL_OUTPUTS - returns model outputs for the passed measuring point and model parameters.
%	GET_RELATIVE_ERROR returns the error and the total error between the passed parameters and the true parameters measured as the difference of the passed and true parameters divided by the true parameters.
%   GET_FACTOR_ERROR returns the error and the total error between the passed parameters and the true parameters measured as the logarithm of the passed parameters divided by the true parameters.
%	GET_ERROR returns the error and the total error between the passed parameters and the true parameters.
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
    
    methods (Static)
        function eta = get_fictitious_measurement_results(model, p, t, var)
        % GET_FICTITIOUS_MEASUREMENT_RESULTS returns fictitious measurement results according to the passed informations.
        %
        % Example:
        %     ETA = MODEL_UTIL.GET_FICTITIOUS_MEASUREMENT_RESULTS(MODEL, P, T, VAR)
        %
        % Input:
        %     MODEL: the model represented as an object whose class
        %            implements the MODEL interface
        %     P: the model parameters according to which the measurement
        %        results are computed
        %     T: the measuring points for which the measurement results are
        %        computed
        %        format: a vector of length n where n is the number of 
        %                measuring points
        %     VAR: the variances of the normal distributed measuring errors
        %          associated with these measuring points
        %          format: a vector of length n where n is the number of 
        %                  measuring points
        %
        % Output:
        %     ETA: the model output for the passed measurements points and
        %          model parameter added to a measurment error according to
        %          the passed variances
        %          format: a vector of length n where n is the number of 
        %                  measuring points
        %
        % see also GET_MODEL_OUTPUTS
        %
        
            out = model_util.get_model_outputs(model, p, t);
            n = length(out);
            eta = zeros(n, 1);

            for i=1:n
                eta(i) =  mvnrnd(out(i), var(i));
            end

        end
        
        function out = get_model_outputs(model, p, t)
        % GET_MODEL_OUTPUTS returns model outputs for the passed measuring point and model parameters.
        %
        % Example:
        %     OUT = MODEL_UTIL.GET_MODEL_OUTPUTS(MODEL, P, T)
        %
        % Input:
        %     MODEL: the model represented as an object whose class
        %            implements the MODEL interface
        %     P: the model parameters according to which the model outputs
        %        are computed
        %     T: the measuring points for which the model outputs are
        %        computed
        %        format: a vector of length n where n is the number of 
        %                measuring points
        %
        % Output:
        %     OUT: the model output for the passed measurements points and
        %          model parameter
        %          format: a vector of length n where n is the number of 
        %                  measuring points
        %
        % see also GET_FICTITIOUS_MEASUREMENT_RESULTS
        %
        
            n = length(t);
            out = zeros(n, 1);

            for i=1:n
                out(i) = model.get_M(p, t(i));
            end
        end
        
        %{
        function [err, total_err] = get_relative_error(p, p_true)
            n = length(p);
            m = length(p_true);
            err = zeros(n, m);
            total_err = zeros(n, 1);
            for i = 1:n
                pi = p(i);
                if iscell(pi)
                    pi = pi{1};
                end
                
                relative_error = (pi - p_true) ./ p_true;
                err(i, 1:m) = relative_error;
                total_err(i) = sum(abs(relative_error));
            end
        end
        
        function [err, total_err] = get_factor_error(p, p_true)
            n = length(p);
            m = length(p_true);
            err = zeros(n, m);
            total_err = zeros(n, 1);
            for i = 1:n
                pi = p(i);
                if iscell(pi)
                    pi = pi{1};
                end
                
                factor_error = log(pi ./ p_true);
                err(i, 1:m) = factor_error;
                total_err(i) = sum(abs(factor_error));
            end
        end
        %}
        
        function [err, total_err] = get_relative_error(p, p_true)  
        % GET_RELATIVE_ERROR returns the error and the total error between the passed parameters and the true parameters measured as the difference of the passed and true parameters divided by the true parameters.
        %
        % Example:
        %     [ERR, TOTAL_ERR] = MODEL_UTIL.GET_RELATIVE_ERROR(P, P_TRUE)
        %
        % Input:
        %     P: a cell array of values which should be compared to P_TRUE
        %        format: a cell array of length n with vectors of length m
        %     P_TRUE: the  value to which the values of P should be
        %             compared to
        %             format: a vectors of length m
        %
        % Output:
        %     ERR: the calculated errors,
        %          ERR(i) = (P(i) - P_TRUE) ./ P_TRUE)
        %          format: a n x m matrix
        %     TOTAL_ERR: the sum of the absolute values of the errors
        %          TOTAL_ERR(i) = abs(sum(ERR(i,:)))
        %          format: a vector of length n
        %
        % see also get_factor_error, get_error
        %
                  
            error_function = @(x, x_true) (x - x_true ./ x_true);
            [err, total_err] = model_util.get_error(p, p_true, error_function);
        end
        
        function [err, total_err] = get_factor_error(p, p_true)
        % GET_FACTOR_ERROR returns the error and the total error between the passed parameters and the true parameters measured as the logarithm of the passed parameters divided by the true parameters.
        %
        % Example:
        %     [ERR, TOTAL_ERR] = MODEL_UTIL.GET_FACTOR_ERROR(P, P_TRUE)
        %
        % Input:
        %     P: a cell array of values which should be compared to P_TRUE
        %        format: a cell array of length n with vectors of length m
        %     P_TRUE: the  value to which the values of P should be
        %             compared to
        %             format: a vectors of length m
        %
        % Output:
        %     ERR: the calculated errors,
        %          ERR(i) = log10(P(i) ./ P_TRUE)
        %          format: a n x m matrix
        %     TOTAL_ERR: the sum of the absolute values of the errors
        %          TOTAL_ERR(i) = abs(sum(ERR(i,:)))
        %          format: a vector of length n
        %
        % see also get_relative_error, get_error
        %
        
            error_function = @(x, x_true) (log10(x ./ x_true));
            [err, total_err] = model_util.get_error(p, p_true, error_function);
        end
        
        function [err, total_err] = get_error(p, p_true, error_function)
        % GET_ERROR returns the error and the total error between the passed parameters and the true parameters.
        %
        % Example:
        %     [ERR, TOTAL_ERR] = MODEL_UTIL.GET_ERROR(P, P_TRUE, ERROR_FUNCTION)
        %
        % Input:
        %     P: a cell array of values which should be compared to P_TRUE
        %        format: a cell array of length n with vectors of length m
        %     P_TRUE: the  value to which the values of P should be
        %             compared to
        %             format: a vectors of length m
        %     ERROR_FUNCTION: the function applied to each comparison to
        %                     quantify the error
        %
        % Output:
        %     ERR: the calculated errors,
        %          ERR(i) = ERROR_FUNCTION(P(i), P_TRUE)
        %          format: a n x m matrix
        %     TOTAL_ERR: the sum of the absolute values of the errors,
        %          TOTAL_ERR(i) = abs(sum(ERR(i,:)))
        %          format: a vector of length n
        %
        % see also get_relative_error, get_factor_error
        %
        
            n = length(p);
            m = length(p_true);
            err = zeros(n, m);
            total_err = zeros(n, 1);
            for i = 1:n
                pi = p(i);
                if iscell(pi)
                    pi = pi{1};
                end
                
                err(i, 1:m) = error_function(pi, p_true);
                total_err(i) = sum(abs(err(i, 1:m)));
            end
        end
    end
    
end

