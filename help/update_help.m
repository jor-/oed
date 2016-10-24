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
clear;
publish('toolbox.m', 'showCode', false);
publish('getting_started.m', 'showCode', false);
publish('product_overview.m', 'showCode', false);
publish('users_guide.m','evalCode', false);
publish('function_overview.m', 'showCode', false);
publish('demo_unidimensional.m');
publish('demo_multidimensional.m');

builddocsearchdb([pwd '/help/html']);

doc('optimal_experimental_design_toolbox');
