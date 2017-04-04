%% Demo: general model
% The purpose of this example is to show how to use an arbitrary MATLAB
% code as a model for the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|>.

%% Create the model object and calculate the needed derivatives by yourself

classdef my_model < model
% This is an example model.
% The needed derivatives are calulated by ourself.
%
% see also MODEL
%

        % ******************* PROPERTIES ******************* %

   properties
      % All class properties are listed here.
      b
      c
   end
   
    methods (Access = public)
        
        % ******************* CONSTRUCTOR ******************* %
        
        function this = my_model(b, c)
        % The setup of the model is done here.
        %
        % Example:
        %     OBJ = MY_MODEL(b, c)
        %
        % Input:
        %     b: an input value for the model setup
        %     c: another input value for the model setup
        %
        % Output:
        %     OBJ: a MY_MODEL object with the passed configurations
        %
        
            this.b = b;
            this.c = c^2;
        end
        
        
        % ******************* MODEL Output ******************* %
        
        function M = get_M(this, p, t)
        % GET_M returns the model output.
        %
        % Example:
        %     M = MY_MODEL.GET_M(P, T)
        %
        % Input:
        %     P: the model parameters
        %     T: the measurement point
        %
        % Output:
        %     M: the model output for the passed parameters P at the measurement point T
        %
        
            M = this.b * p(1)^2 * t(1) + this.c * p(2) * t(2);
        end
        
        function dp_M = get_dp_M(this, p, t)
        % GET_DP_M returns the derivative of the model output with respect to the model parameters.
        %
        % Example:
        %     M = MY_MODEL.GET_DP_M(P, T)
        %
        % Input:
        %     P: the model parameters
        %     T: the measurement point
        %
        % Output:
        %     M: the derivative of the model output with respect to the model parameters
        %        for the passed parameters P at the measurement T
        %
        
            dp_M = [this.b * 2 * p(1) * t(1); this.c * t(2)];
        end
        
        function dpdp_M = get_dpdp_M(this, p, t)
        % GET_DPDP_M returns the second derivative of the model output with respect to the model parameters.
        %
        % Example:
        %     M = MY_MODEL.GET_DPDP_M(P, T)
        %
        % Input:
        %     P: the model parameters
        %     T: the measurement point
        %
        % Output:
        %     M: the second derivative of the model output with respect to the model parameters
        %        for the passed parameters P at the measurement T
        %
        
            dpdp_M = [this.b * 2 * t(1), 0; 0, 0];
        end
        
    end
    
end

b = 1;
c = 2;
model_object = my_model(b, c);
