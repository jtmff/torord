% This script compares the behaviour of fixed-Cl and dynamic-Cl versions of
% ToR-ORd after 1k beats (I.e., this does not really showcase the
% instability of the fixed-Cl model, it serves mainly to show general
% similarity of action potential and calcium transient at the time when the
% model is usually evaluated. 
%% Setting parameters
clear 

param.bcl = 1000; % basic cycle length in ms
param.model = @model_Torord_dynCl; 


params(1:2) = param;
params(2).model = @model_Torord;


options = []; % parameters for ode15s - usually empty
beats = 1000; % number of beats
ignoreFirst = 999; %999;% beats - 1; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

X0{1} = getStartingState('Torord_endo_dynCl'); 
X0{2} = getStartingState('Torord_endo'); 



%% Simulation and extraction of outputs

parfor i = 1:2
    [time{i}, X{i}] = modelRunner(X0{i}, options, params(i), beats, ignoreFirst);
    currents{i} = getCurrentsStructure(time{i}, X{i}, params(i), 0);
end


%% Plotting membrane potential and calcium transient, as well as chloride

figure(1);
plot(currents{1}.time, currents{1}.V, currents{2}.time, currents{2}.V,':', 'LineWidth', 1.5, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
xlim([0 500]);
legend({'ToR-ORd-dynCl','ToR-ORd'}, 'Location','east');

figure(2);
plot(currents{1}.time, currents{1}.Cai, currents{2}.time, currents{2}.Cai,':', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Ca_i (mmol/L)');
legend({'ToR-ORd-dynCl','ToR-ORd'}, 'Location','east');



