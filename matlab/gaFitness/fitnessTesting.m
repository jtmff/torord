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

function [errors] = fitnessTesting(logMultipliers)
%FITNESSTESTING Returns fitness of a ToR-ORd model with modifications
% The fitness in variable errors is a 2-element vector, with errors(1) coding for squared error of
% APD versus 400, and errors(2) codes for squared error of CaT versus 200
% nM (2e-4 mM).

%% Making log-multipliers multipliers again.
currentMultipliers = exp(logMultipliers);

%% Setting parameters (multipliers are used here
param.bcl = 1000;
param.cellType = 0; % endocardial
param.ICaL_Multiplier = currentMultipliers(1);
param.IKr_Multiplier = currentMultipliers(2);
param.IKs_Multiplier = currentMultipliers(3);
param.INaL_Multiplier = currentMultipliers(4);
param.Jup_Multiplier = currentMultipliers(5);

% The next parameter setting is commented out at the moment - if uncommented, it slows the runtime by ca.
% 5-10%, but it makes sure no beat is simulated for much more than 5 seconds. This
% is key when the GA can produce nonsensical models (e.g. when currents can
% go crazy high) that crash - such crashes take ages to happen when ode15s
% is used, and they are extremely costly, as they eat a processor core for
% a long time and they can stall a generation of the GA too. In our
% experience, this is not a problem in the sample GA given here, but it did
% occur to us when we were making the ToR-ORd model.
% param.maxTimePerBeat = 0.5; 

options = [];
%% Standard setting of parameter simulations and running the simulation

beats = 50; % number of beats - very short for demonstration purposes
ignoreFirst = beats - 1; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

X0 = getStartingState('Torord_endo');

[time, X, parameters] = modelRunner(X0, options, param, beats, ignoreFirst);
currents = getCurrentsStructure(time, X, param, 0);

if parameters.isFailed == 0 % If this condition holds, the simulation finished fine.   
    %% Computing biomarkers and returning fitness based on them.
    apd90 = DataReporter.getAPD(currents.time, currents.V, 0.9);
    CaTamplitude = max(currents.Cai) - min(currents.Cai);
    
    errors(1) = (apd90 - 400)^2;
    errors(2) = (CaTamplitude - 2e-4)^2; % note that even though the magnitude is completely different from the 1st dimension, this is not a problem, as the 2 dimensions are treated separately.
else % If simulation terminated prematurely, it's most likely useless
    errors(1) = Inf;
    errors(2) = Inf;
end
end

