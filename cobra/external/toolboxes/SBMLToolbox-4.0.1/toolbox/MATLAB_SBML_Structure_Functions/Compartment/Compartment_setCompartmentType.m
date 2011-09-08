function SBMLCompartment = Compartment_setCompartmentType(SBMLCompartment, compartmentType)
% SBMLCompartment = Compartment_setCompartmentType(SBMLCompartment, compartmentType)
%
% Takes
%
% 1. SBMLCompartment, an SBML Compartment structure
% 2. compartmentType; a string representing the compartmentType to be set
%
% Returns
%
% 1. the SBML Compartment structure with the new value for the compartmentType attribute
%

%<!---------------------------------------------------------------------------
% This file is part of SBMLToolbox.  Please visit http://sbml.org for more
% information about SBML, and the latest version of SBMLToolbox.
%
% Copyright (C) 2009-2011 jointly by the following organizations: 
%     1. California Institute of Technology, Pasadena, CA, USA
%     2. EMBL European Bioinformatics Institute (EBML-EBI), Hinxton, UK
%
% Copyright (C) 2006-2008 jointly by the following organizations: 
%     1. California Institute of Technology, Pasadena, CA, USA
%     2. University of Hertfordshire, Hatfield, UK
%
% Copyright (C) 2003-2005 jointly by the following organizations: 
%     1. California Institute of Technology, Pasadena, CA, USA 
%     2. Japan Science and Technology Agency, Japan
%     3. University of Hertfordshire, Hatfield, UK
%
% SBMLToolbox is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation.  A copy of the license agreement is provided
% in the file named "LICENSE.txt" included with this software distribution.
%----------------------------------------------------------------------- -->
























%get level and version and check the input arguments are appropriate

[level, version] = GetLevelVersion(SBMLCompartment);

if isfield(SBMLCompartment, 'compartmentType')
	if ~ischar(compartmentType)
		error('compartmentType must be character array') ;
	else
		SBMLCompartment.compartmentType = compartmentType;
	end;
else
	error('compartmentType not an attribute on SBML L%dV%d Compartment', level, version);
end;

