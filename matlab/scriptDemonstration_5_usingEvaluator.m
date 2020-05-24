%% This script shows how to run an automated model evaluator. In this case, reportParams.reportType = 'full' determines that all Matlab-based criteria are computed and plotted
% the evaluator runs the basic version of the model. To run the dynCl
% version, please uncomment lines 10 and 18.
clear

% First, we define which model is to be evaluated, using the same type of
% parameter structure as before.
param.bcl = 1000;
param.model = @model_Torord;
% param.model = @model_Torord_dynCl; % Uncomment for dynCl version 

folderOut = 'data';  % output folder
fnameSave = 'modelReport'; % filename with results of simulations ran using the given model
reportParams.reportType = 'full'; % Plotting both calibration and validation. Otherwise, this can be 'calibration', 'validation'.

switches.parameter = param; % This way, we tell the evaluator which model to use.
switches.X0 = getStartingState('Torord_endo'); % And this specifies the starting state
% switches.X0 = getStartingState('Torord_endo_dynCl'); % Uncomment for dynCl version

%% Running evaluation simulations (takes a while, can be around an hour or so)
DataReporter.runEvaluationSimulations( switches, folderOut, fnameSave, reportParams )

%% Generating reports and plots.
fileSimulation = [folderOut '/' fnameSave];
fileExperiment = 'data/oliTraces'; % Path to file containing experimental data on human AP morphology (used to show/assess agreement with model's AP trace)
folderOutReport = 'data/modelReport_full'; % Path to the folder with produced report.
DataReporter.makeReport(fileSimulation, fileExperiment, folderOutReport, reportParams)
