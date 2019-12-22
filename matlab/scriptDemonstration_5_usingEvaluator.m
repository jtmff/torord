%     Cardiac model ToR-ORd
%     Copyright (C) 2019 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


%% This script shows how to run an automated model evaluator. In this case, reportParams.reportType = 'full' determines that all Matlab-based criteria are computed and plotted
clear

% First, we define which model is to be evaluated, using the same type of
% parameter structure as before.
param.bcl = 1000;
param.model = @model_Torord;

folderOut = 'data';  % output folder
fnameSave = 'm12Report'; % filename with results of simulations ran using the given model
reportParams.reportType = 'full'; % Plotting both calibration and validation. Otherwise, this can be 'calibration', 'validation'.

switches.parameter = param; % This way, we tell the evaluator which model to use.
switches.X0 = getStartingState('Torord_endo'); % And this specifies the starting state

%% Running evaluation simulations (takes a while, can be around an hour or so)
DataReporter.runEvaluationSimulations( switches, folderOut, fnameSave, reportParams )

%% Generating reports and plots.
fileSimulation = [folderOut '/' fnameSave];
fileExperiment = 'data/oliTraces'; % Path to file containing experimental data on human AP morphology (used to show/assess agreement with model's AP trace)
folderOutReport = 'data/m12Report_full'; % Path to the folder with produced report.
DataReporter.makeReport(fileSimulation, fileExperiment, folderOutReport, reportParams)
