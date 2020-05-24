%% This gives an example of a multicriterial genetic algorithm which refits the ToR-ORd model to aim for 
% APD 400 (greater than in control) and CaT amplitude of 2e-4 (200 nM),
% which is smaller than in control. This is of course just an artificial
% example with no particular meaning. If you're unfamiliar with this method
% in the context of cardiac models, feel free to have a look at our blog at
% TODO link when published.
%
% The parameters that are evolved are conductances of: ICaL, IKs, IKr,
% INaL, and Jup multiplier (SERCA reuptake). I.e., we have 5 input
% parameters, and a 2-element fitness.

clear

addpath gaFitness; % this folder contains fitness functions

disp('starting fitting');

% Parameters for the optimizer:
nGenerations = 10; % number of optimization generations (cycles of crossover, mutation, and selection)
timeLimit = 60 * 60; % rough time limit of the optimization in seconds. However, the generation that happens when the limit is reached is finished.
popSize = 50; % Number of creatures in the population.

% The runtime of the optimization is roughly proportional to X * popSize * nGenerations,
% where X is the runtime of a single fitness (worth measuring this using
% tic-toc). The runtime is then cut down by the number of cores used. The
% scaling of this is good and performance gain is substantial, as long as
% the fitness computation takes a nontrivial amount of time (>1s), which it
% does almost always with cardiac models.

% Here are options set - in addition to the ones above, 'PlotFcn' selects
% type of progress visualisation (plotting Pareto front in this case), and
% 'UseParallel' says that all available cores are to be used. If you run
% parpool(X) before this line, it specifies how many cores are
% to be used (e.g. less than maximum if the computer is to be used for
% other work).
options = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto, 'UseParallel',1, 'Generations', nGenerations, 'TimeLimit', timeLimit, 'PopulationSize', popSize);

% We can restrict the multipliers of curents. Here, we restrict them to 50%-200% of the current. As we're using
% log-multipliers (see Appendix/SupMat 2 for rationale), these correspond
% to log(0.5) ~ -0.7 and log(2) ~ 0.7.
lowerBound = -0.7 * ones(5,1);  
upperBound = 0.7 * ones(5,1);  

% The multicriterial GA is started (this takes a while even for the small
% population and short fitness runtime. It also does not really converge -
% more creatures and generations would be needed.

% @fitnessTesting is a handle to the fitness function.
[coeffEstimate,FVAL,EXITFLAG,OUTPUT,POPULATION,SCORE] = gamultiobj(@fitnessTesting, 5,[],[],[],[], lowerBound, upperBound, options);

% coeffEstimate is encoding of the Pareto front, where each line
% corresponds to a creature and each column
save data/gaOutputs.mat; 

%% Storing currents and errors of points on the Pareto front
load data/gaOutputs.mat;
for iCreature = 1:size(coeffEstimate, 1)
    [ errorvalues{iCreature}, currentsModels{iCreature}] = fitnessTestingCurrents(coeffEstimate(iCreature,:)); % fitnessTestingCurrents is copy-pasted fitnessTesting, which furthermore returns the structure of currents, so we can look at how the optimized cell looks.
end

%% Plotting outputs of the models. The models best in one criterion are usually quite bad in the other one, so it's mainly about whether the models that are pretty good in both criteria are good enough.
for iCreature = 1:size(coeffEstimate, 1)
    
    apd90 = DataReporter.getAPD(currentsModels{iCreature}.time, currentsModels{iCreature}.V, 0.9);
    CaTamplitude = max(currentsModels{iCreature}.Cai) - min(currentsModels{iCreature}.Cai);
    
    figure(iCreature); clf
    
    subplot(2,1,1);
    plot(currentsModels{iCreature}.time, currentsModels{iCreature}.V);
    xlabel('Time (ms)');
    ylabel('Mem. pot. (mV)');
    xlim([0 600]);
    title(['APD90: ' num2str(apd90) ' ms']);
    
    subplot(2,1,2);
    plot(currentsModels{iCreature}.time, currentsModels{iCreature}.Cai);
    xlabel('Time (ms)');
    ylabel('Cai (mM)');
    xlim([0 600]);
    title(['CaT amplitude: ' num2str(CaTamplitude*10^6) ' nM']);
end

%% COMMENT: When large population is computed with long generations, it's often worth storing each generation (in case of power outage, or researcher wanting
% to look at non-final generations, to see how the algorithm converges). To do this, outputFcn (set in gaoptimset) can be defined to store population 
% at the end of each generation.
%
% That may look e.g. like:
%
% function [state,options,optchanged] = saveGeneration(options,state,flag)
% generation = state.Generation;
% population = state.Population;
% score = state.Score;
% optchanged = false;
% save(['intermediateGAstates/' num2str(generation)], 'population','score');
% end
%
% % Note that the folder intermediateGAstates must exist first.
