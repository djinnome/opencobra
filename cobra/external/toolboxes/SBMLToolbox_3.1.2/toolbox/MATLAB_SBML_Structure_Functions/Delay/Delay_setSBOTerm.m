function SBMLDelay = Delay_setSBOTerm(SBMLDelay, sboTerm)
%
%   Delay_setSBOTerm 
%             takes  1) an SBMLDelay structure 
%             and    2) an integer representing the sboTerm to be set
%
%             and returns 
%               the Delay with the sboTerm set
%
%       SBMLDelay = Delay_setSBOTerm(SBMLDelay, sboTerm)

%  Filename    :   Delay_setSBOTerm.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: Delay_setSBOTerm.m 7155 2008-06-26 20:24:00Z mhucka $
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
if (~isstruct(SBMLDelay))
    error(sprintf('%s\n%s', ...
      'Delay_setSBOTerm(SBMLDelay, sboTerm)', ...
      'argument must be an SBML Delay structure'));
end;
 
[sbmlLevel, sbmlVersion] = GetLevelVersion(SBMLDelay);

if (~isSBML_Delay(SBMLDelay, sbmlLevel, sbmlVersion))
  error(sprintf('%s\n%s', ...
    'Delay_setSBOTerm(SBMLDelay, sboTerm)', ...
    'first argument must be an SBML Delay structure'));
elseif (~isIntegralNumber(sboTerm))
    error(sprintf('%s\n%s', ...
      'Delay_setSBOTerm(SBMLDelay, sboTerm)', ...
      'second argument must be an integer representing the sboTerm'));
elseif (sbmlLevel ~= 2 || sbmlVersion ~= 3)
    error(sprintf('%s\n%s', ...
      'Delay_setSBOTerm(SBMLDelay, sboTerm)',  ...
      'sboTerm field only in level 2 version 3 model'));    
end;

SBMLDelay.sboTerm = sboTerm;
