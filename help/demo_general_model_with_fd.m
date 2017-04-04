%% Demo: general model with finite differences
% The purpose of this example is to show how to use an arbitrary MATLAB
% code as a model for the <matlab:doc('optimal_experimental_design_toolbox')
% |Optimal Experimental Design Toolbox|>.

%% Create the model object and calculate the needed derivatives by finite differences

classdef my_model < model_fd
% This is an example model.
% The needed derivatives are calulated using finite differences.
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
        
    end
    
end

b = 1;
c = 2;
model_object = my_model(b, c);
