% Optimal Experimental Design Toolbox
% Version 1.1.0 07-Apr-2017
%
% Solver.
%   solver             - allows to calculate and optimize the quality of experimental designs and allows to calculate a parameter estimation resulting from accomplished measurements.
%   solver_options     - represents the options for a solver object.
%   solver_edo_options - represents the options for the solver of the experimental design optimization problem.
%   solver_po_options  - represents the options for the solver of the parameter optimization problem.
%
% Criteria.
%   criterion          - represents an interface for a quality criterion.
%   criterion_A        - is the quality criterion that uses the average variance as the quality.
%
% Models.
%   model              - represents an interface for a model.
%   model_explicit     - implements the model interface and provides the function value and the first and second derivatives with respect to the parameters of an explicitly given (symbolic) model function.
%   model_ivp          - implements the model interface and provides the solution of an initial value problem and its first and second derivatives with respect to the parameters.
%   model_fd           - implements the model interface and provides the first and second derivatives with respect to the parameters by finite differences approximations.
%   model_composed     - implements the model interface and provides the function value and the first and second derivatives with respect to the parameters of an explicitly given (symbolic) model function including an inner model.
%
% Utilities.
%   util               - represents an utility class with static auxiliary methods for the toolbox.
%
% Application examples.
%   model_C2           - models the concentration C of the suspended sediment above the marsh surface with two parameters.
%   model_C3           - models the concentration C of the suspended sediment above the marsh surface with three parameters.
%   model_util         - represents an utility class with static auxiliary methods for working with model classes.

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
