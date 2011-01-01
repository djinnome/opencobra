function reaction = Reaction_addReactant(SBMLReaction, SBMLReactant)
%
%   Reaction_addReactant 
%             takes  1) an SBMLReaction structure 
%             and    2) an SBMLReactant structure
%
%             and returns 
%               the reaction with the reactant added
%
%       reaction = Reaction_addReactant(SBMLReaction, SBMLReactant)

%  Filename    :   Reaction_addReactant.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: Reaction_addReactant.m 7155 2008-06-26 20:24:00Z mhucka $
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
if (~isstruct(SBMLReaction))
  error(sprintf('%s\n%s', ...
    'Reaction_addReactant(SBMLReaction, SBMLReactant)', ...
    'first argument must be an SBML Reaction structure'));
end;
 
[sbmlLevel, sbmlVersion] = GetLevelVersion(SBMLReaction);

if (~isSBML_Reaction(SBMLReaction, sbmlLevel, sbmlVersion))
    error(sprintf('%s\n%s', 'Reaction_addReactant(SBMLReaction, SBMLReactant)', 'first argument must be an SBML reaction structure'));
elseif (~isSBML_SpeciesReference(SBMLReactant, sbmlLevel, sbmlVersion))
    error(sprintf('%s\n%s\n of the same SBML level, namely level %u', 'Reaction_addReactant(SBMLReaction, SBMLReactant)', 'second argument must be an SBML reactant structure', sbmlLevel));
end;

numberReactants = length(SBMLReaction.reactant);

SBMLReaction.reactant(numberReactants+1) = SBMLReactant;

reaction = SBMLReaction;

