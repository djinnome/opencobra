function units = ParameterRule_getUnits(SBMLParameterRule)
%
%   ParameterRule_getUnits 
%             takes an SBMLParameterRule structure 
%
%             and returns 
%               the units of the parameterRule as a string
%
%       units = ParameterRule_getUnits(SBMLParameterRule)

%  Filename    :   ParameterRule_getUnits.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: ParameterRule_getUnits.m 7155 2008-06-26 20:24:00Z mhucka $
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



% check that input is correct
if (~isstruct(SBMLParameterRule))
  error(sprintf('%s\n%s', ...
    'ParameterRule_getUnits(SBMLParameterRule)', ...
    'first argument must be an SBML ParameterRule structure'));
end;
 
[sbmlLevel, sbmlVersion] = GetLevelVersion(SBMLParameterRule);

if (~isSBML_ParameterRule(SBMLParameterRule, sbmlLevel, sbmlVersion))
    error(sprintf('%s\n%s', 'ParameterRule_getUnits(SBMLParameterRule)', 'argument must be an SBML parameterRule structure'));
end;

units = SBMLParameterRule.units;
