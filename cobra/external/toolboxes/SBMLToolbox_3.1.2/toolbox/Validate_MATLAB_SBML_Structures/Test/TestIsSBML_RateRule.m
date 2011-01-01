function fail = TestIsSBML_RateRule

%  Filename    :   TestIsSBML_RateRule.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: TestIsSBML_RateRule.m 7155 2008-06-26 20:24:00Z mhucka $
%  $Source v $
%
%<!---------------------------------------------------------------------------
% This file is part of SBMLToolbox.  Please visit http://sbml.org for more
% information about SBML, and the latest version of SBMLToolbox.
%
% Copyright 2005-2007 California Institute of Technology.
% Copyright 2002-2005 California Institute of Technology and
%                     Japan Science and Technology Corporation.
% 
% This library is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation.  A copy of the license agreement is provided
% in the file named "LICENSE.txt" included with this software distribution.
% and also available online as http://sbml.org/software/sbmltoolbox/license.html
%----------------------------------------------------------------------- -->


rr_l1 = struct('typecode', {'SBML_PARAMETER_RULE'}, 'notes', {''}, 'annotation', {''}, 'type', ...
    {'rate'}, 'formula', {''}, 'variable', {''}, 'species', {''}, 'compartment', {''}, 'name', {''}, 'units', {''});

rr_l2 = struct('typecode', {'SBML_RATE_RULE'}, 'metaid', {''}, 'notes', {''}, 'annotation', {''},  ...
    'formula', {''}, 'variable', {''}, 'species', {''}, 'compartment', {''}, 'name', {''}, 'units', {''});

rr_l2v2 = struct('typecode', {'SBML_RATE_RULE'}, 'metaid', {''}, 'notes', {''}, 'annotation', {''}, 'sboTerm', {''}, ...
    'formula', {''}, 'variable', {''}, 'species', {''}, 'compartment', {''}, 'name', {''}, 'units', {''});

fail = TestFunction('isSBML_RateRule', 2, 1, rr_l1, 1, 1);
fail = fail + TestFunction('isSBML_RateRule', 3, 1, rr_l1, 1, 1, 1);
fail = fail + TestFunction('isSBML_RateRule', 3, 1, rr_l1, 1, 2, 1);
fail = fail + TestFunction('isSBML_RateRule', 2, 1, rr_l2, 2, 1);
fail = fail + TestFunction('isSBML_RateRule', 3, 1, rr_l2, 2, 1, 1);
fail = fail + TestFunction('isSBML_RateRule', 3, 1, rr_l2v2, 2, 2, 1);










