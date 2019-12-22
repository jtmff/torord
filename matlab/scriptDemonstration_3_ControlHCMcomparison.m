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


%% A script showing how two different models can be easily compared. Here,
% we plot control ToR-ORd and the HCM variant in control condition and with
% 50% IKr block
% Slow pacing (0.25 Hz) is used to facilitate EADs. 

% Setting parameters
clear 
% param is the default model parametrization here
param.bcl = 4000;
param.model = @model_Torord;
param.IKr_Multiplier = 1; 
param.verbose = 1; % It is possible to use verbose output even when parfor is used for simulations (the numbers from threads are mixed together, but they still given an idea of the progress).

params(1:4) = param; % first 2 are control, last 2 HCM. Every second model in these conditions has 50% IKr block. By passing a different model function handle as the 'model' parameter, we can easily simulate different model files within the same framework.
params(3).model = @model_Torord_HCM; params(4).model = @model_Torord_HCM;
params(2).IKr_Multiplier = 0.5; params(4).IKr_Multiplier = 0.5;


options = [];
beats = 1000;
ignoreFirst = beats - 1;

%% Simulation and output extraction

parfor i = 1:length(params)
    X0 = getStartingState('Torord_endo'); 
    [time{i}, X{i}] = modelRunner(X0, options, params(i), beats, ignoreFirst);
    currents{i} = getCurrentsStructure(time{i}, X{i}, params(i), 0);
end

%% Plotting APs
figure(11); clf
for i = 1:length(params)
    hold on
    plot(currents{i}.time, currents{i}.V);
    hold off
end

xlim([0 1500]);

legend('ToR-ORd', 'ToR-ORd + 50% I_{Kr} block','ToR-ORd HCM', 'ToR-ORd HCM + 50% I_{Kr} block');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
% Please note that while the IKr-blocked models are the identical to Figure
% 8E, the other two traces are not the same as in Figure 8C, as that (as
% well as 8D) were paced at 1000 ms bcl.

