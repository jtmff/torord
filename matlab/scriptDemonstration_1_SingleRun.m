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



%% This is a simple script which runs the control endocardial model for 100
% beats and plots membrane potential and calcium transient.

% Setting parameters
clear 

% Param is the structure of model parameters that the user may wish to
% change compared to default simulation. The full list is given in the
% function ORdRunner, and it mainly includes cell type, current
% multipliers, extracellular ionic concentrations, or fraction of NCX and ICaL
% localisation in junctional subspace.
param.bcl = 1000; % basic cycle length in ms
param.model = @model_Torord; % which model is to be used - right now, use model12. In general, any model with the same format of inputs/outputs as model12 may be simulated, which is useful when current formulations are changed within the model code, etc.
param.verbose = true; % printing numbers of beats simulated.

options = []; % parameters for ode15s - usually empty
beats = 100; % number of beats
ignoreFirst = beats - 1; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

X0 = getStartingState('Torord_endo'); % starting state - can be also m12_mid or m12_epi for midmyocardial or epicardial cells respectively.

%% Simulation and extraction of outputs

% The structure param and other variables are passed to ORdRunner, which is
% an interface between user and the simulation code itself (which is in
% model12.m). The ORdRunner unpacks the structure of parameters given by
% the users, sets undefined parameters to default, and sends all that to
% @model12.

% time, X are cell arrays corresponding to stored beats (if 1 beat is
% simulated, this is 1-by-1 cell still), giving time vectors and state
% variable values at corresponding time points.
[time, X] = modelRunner(X0, options, param, beats, ignoreFirst);

% A structure of currents is computed from the state variables (see the
% function code for a list of properties extracted - also, hitting Tab
% following typing 'currents.' lists all the fields of the structure). Some
% state variables are also stored in a named way (time, V, Cai, Cass) so
% that the user can do most of necessary plotting simply via accessing the
% structure currents as shown below. 
currents = getCurrentsStructure(time, X, param, 0);


%% Plotting membrane potential and calcium transient
figure(1);
plot(currents.time, currents.V);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
xlim([0 500]);

figure(2);
plot(currents.time, currents.Cai);
xlabel('Time (ms)');
ylabel('Ca_i (mM)');
xlim([0 500]);