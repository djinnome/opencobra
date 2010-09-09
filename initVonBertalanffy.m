%setup the paths to the data, scripts and functions, load data, +/- generate docs

%clear workspace
clear

%find the absolute path of this COBRA toolbox extension
global vonBdir
vonBdir=which('generateVonBertalanffyDocumentation');
if isempty(vonBdir)
    error('Current directory must be */vonBertalanffy')
end
vonBdir=vonBdir(1:end-length('/generateVonBertalanffyDocumentation.m'));

load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'alberty2006' filesep 'Alberty2006.mat']);

load([vonBdir filesep 'setupThermoModel' filesep  'experimentalData' filesep 'alberty2006' filesep 'metAbbrAlbertyAbbr.mat']);

%group contrbuion data for E. coli
load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'groupContribution' filesep 'metGroupCont_Ecoli_iAF1260.mat']);

if isempty(which('pHbalanceProtons'))
    fprintf('Add the toolbox to the matlab path.\n')
    %add the toolbox path to the matlab path
    addpath(genpath(vonBdir));
    savepath;
end

if ~exist([vonBdir filesep 'doc'])
    fprintf('Generating html documentation... \n')
    generateVonBertalanffyDocumentation;
    fprintf('...documentation complete.\n')
end
