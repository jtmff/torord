% A script demonstrating how updates to the ionic activity coefficients involving decadic rather natural logarithm (see our comment/annotation at the
% Elife article https://elifesciences.org/articles/48890) maintain a very similar-to-the-original
% behaviour both for the fixed-Cl and dynamic-Cl variant of the model.
%% Setting parameters
clear 
addpath otherModelCodes;  % This is where the "updatedPhiCaL" model codes are placed for params(2) and params(4). Briefly, the ionic activity coefficients are computed using 10^ rather than exp(), and conductances are scaled accordingly. See the comment attached to our article for full details: https://elifesciences.org/articles/48890.

% param is the default model parametrization here
param.bcl = 1000;
param.model = @model_Torord;
param.verbose = 1; % It is possible to use verbose output even when parfor is used for simulations (the numbers from threads are mixed together, but they still given an idea of the progress).

params(1:4) = param; % first 2 are control, last 2 HCM. Every second model in these conditions has 50% IKr block. By passing a different model function handle as the 'model' parameter, we can easily simulate different model files within the same framework.
params(2).model = @model_Torord_updatedPhical;
params(3).model = @model_Torord_dynCl;
params(4).model = @model_Torord_dynCl_updatedPhiCaL;


options = [];
beats = 1000;
ignoreFirst = beats - 1;

%% Simulation and output extraction

parfor i = 1:length(params)
    if (i <= 2)
        X0 = getStartingState('Torord_endo'); 
    else
        X0 = getStartingState('Torord_endo_dynCl'); 
    end
    [time{i}, X{i}] = modelRunner(X0, options, params(i), beats, ignoreFirst);
    currents{i} = getCurrentsStructure(time{i}, X{i}, params(i), 0);
end

%% Plotting AP, CaT, and ICaL, showing how the updated PhiCaL versions are near-identical to the original models (with the dynamic Cl model showing only very minor difference to the fixed Cl model).
figure(1); clf
for i = 1:length(params)
    hold on
    plot(currents{i}.time, currents{i}.V);
    hold off
end

xlim([0 500]);

legend('ToR-ORd', 'ToR-ORd + updated PhiCaL','ToR-ORd dynCl', 'ToR-ORd dynCl + updated PhiCaL');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');

figure(2); clf
for i = 1:length(params)
    hold on
    plot(currents{i}.time, currents{i}.Cai);
    hold off
end

xlim([0 500]);

legend('ToR-ORd', 'ToR-ORd + updated PhiCaL','ToR-ORd dynCl', 'ToR-ORd dynCl + updated PhiCaL');
xlabel('Time (ms)');
ylabel('Cai');

figure(3); clf
for i = 1:length(params)
    hold on
    plot(currents{i}.time, currents{i}.ICaL);
    hold off
end

xlim([0 500]);

legend('ToR-ORd', 'ToR-ORd + updated PhiCaL','ToR-ORd dynCl', 'ToR-ORd dynCl + updated PhiCaL');
xlabel('Time (ms)');
ylabel('ICaL');
