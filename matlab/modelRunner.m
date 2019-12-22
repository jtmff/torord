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
function [ time, X, parameters ] = modelRunner( X0, options, parameters, beats, ignoreFirst)
% A function for simulating models with various parameters. It serves as
% the interface between user's scripts (defining parameters) and simulation
% core. Here, the structure of parameters is unpacked and passed to the model
% INPUTS:
% X0 - starting state.
%
% options - options to ode15s (may be left empty if not sure what to do with it).
%
% parameters - passed to the model and specify, e.g. scaling of conductances or which implementation of a current is used.
%
% beats - the number of beats in the simulation
%
% ignoreFirst - specifies how many action potentials are to be ignored out of beats. This serves to save memory
%
% OUTPUTS:
% time - a cell array; the i-th element gives the timeline of the i-th
% action potential
%
% X - a cell array; the i-th element gives the matrix of state variables (#rows = length of
% time vector, columns are state variables).
%
% parameters - the same parameters that were passed to this function

%% Parameters are set here
% defaults which may be overwritten

% Cell type (endo/epi/mid)
cellType = 0; if (isfield(parameters,'cellType')) cellType = parameters.cellType;end

% If the number of simulated beat is to be printed out.
verbose = false; if (isfield(parameters,'verbose')) verbose = parameters.verbose; end

% Extracellular concentrations
nao = 140; if (isfield(parameters,'nao')) nao = parameters.nao; end
cao = 1.8;if (isfield(parameters,'cao')) cao = parameters.cao; end
ko = 5;if (isfield(parameters,'ko')) ko = parameters.ko; end

% Localization of ICaL and NCX: the fraction in junctional subspace
ICaL_fractionSS = 0.8; if (isfield(parameters, 'ICaL_fractionSS')) ICaL_fractionSS = parameters.ICaL_fractionSS; end
INaCa_fractionSS = 0.35; if (isfield(parameters, 'INaCa_fractionSS')) INaCa_fractionSS = parameters.INaCa_fractionSS; end

% Current multipliers
INa_Multiplier = 1; if (isfield(parameters,'INa_Multiplier')) INa_Multiplier = parameters.INa_Multiplier; end
ICaL_Multiplier = 1; if (isfield(parameters,'ICaL_Multiplier')) ICaL_Multiplier = parameters.ICaL_Multiplier; end
Ito_Multiplier = 1; if (isfield(parameters,'Ito_Multiplier')) Ito_Multiplier = parameters.Ito_Multiplier; end
INaL_Multiplier = 1; if (isfield(parameters,'INaL_Multiplier')) INaL_Multiplier = parameters.INaL_Multiplier; end
IKr_Multiplier = 1; if (isfield(parameters,'IKr_Multiplier')) IKr_Multiplier = parameters.IKr_Multiplier; end
IKs_Multiplier = 1; if (isfield(parameters,'IKs_Multiplier')) IKs_Multiplier = parameters.IKs_Multiplier; end
IK1_Multiplier = 1; if (isfield(parameters,'IK1_Multiplier')) IK1_Multiplier = parameters.IK1_Multiplier; end
IKb_Multiplier = 1; if (isfield(parameters,'IKb_Multiplier')) IKb_Multiplier = parameters.IKb_Multiplier; end
INaCa_Multiplier = 1; if (isfield(parameters,'INaCa_Multiplier')) INaCa_Multiplier = parameters.INaCa_Multiplier; end
INaK_Multiplier = 1; if (isfield(parameters,'INaK_Multiplier')) INaK_Multiplier = parameters.INaK_Multiplier; end
INab_Multiplier = 1; if (isfield(parameters,'INab_Multiplier')) INab_Multiplier = parameters.INab_Multiplier; end
ICab_Multiplier = 1; if (isfield(parameters,'ICab_Multiplier')) ICab_Multiplier = parameters.ICab_Multiplier; end
IpCa_Multiplier = 1; if (isfield(parameters,'IpCa_Multiplier')) IpCa_Multiplier = parameters.IpCa_Multiplier; end
ICaCl_Multiplier = 1; if (isfield( parameters, 'ICaCl_Multiplier')) ICaCl_Multiplier =  parameters.ICaCl_Multiplier; end
IClb_Multiplier = 1; if (isfield( parameters, 'IClb_Multiplier')) IClb_Multiplier =  parameters.IClb_Multiplier; end
Jrel_Multiplier = 1; if (isfield(parameters,'Jrel_Multiplier')) Jrel_Multiplier = parameters.Jrel_Multiplier; end
Jup_Multiplier = 1; if (isfield(parameters,'Jup_Multiplier')) Jup_Multiplier = parameters.Jup_Multiplier; end

% An array of extra parameters (if user-defined, this can be a cell array
% as well), passing other parameters not defined otherwise in this script
extraParams = []; if (isfield(parameters,'extraParams')) extraParams = parameters.extraParams; end

% There may be parameters defining clamp behaviour
vcParameters = []; if (isfield(parameters, 'vcParameters')) vcParameters = parameters.vcParameters; end
apClamp = []; if (isfield(parameters,'apClamp')) apClamp = parameters.apClamp; end

% Stimulus parameters
stimAmp = -53; if (isfield(parameters,'stimAmp')) stimAmp = parameters.stimAmp; end
stimDur = 1; if (isfield(parameters,'stimDur')) stimDur = parameters.stimDur; end

% The parameter maxTimePerBeat determines the maximum runtime of one beat
% simulation - and kills the computation prematurely if it's not the case. It slows
% the runtime by ca. 10%, but it makes sure no simulation takes ages (which
% can happen when a genetic algorithm/population of models produces a
% nonsensical model which is unstable and crashes, but ode15s keeps making
% dt smaller forever, so it can take 4 hours before the crash finally
% occurs). This may be important when dysfunctional models stall a
% generation in an optimizer.
maxTimePerBeat = Inf;
if (isfield(parameters,'maxTimePerBeat'))  % Maximum time of simulation per beat can be specified - only to be used when encountering unstable models that would crash anyway, but would take hours to do that.
    maxTimePerBeat = parameters.maxTimePerBeat;
    
    xoverFcn = @(T, Y, varargin) myEventFunction(T, Y, maxTimePerBeat, varargin ); % Many thanks to https://uk.mathworks.com/matlabcentral/answers/36683-timeout-using-parfor-loop-and-ode15s#answer_130923 for suggesting this trick, as well as to https://uk.mathworks.com/matlabcentral/profile/authors/6063611-daniel-m for pointing me to it.
    
    options = odeset(options, 'Events', xoverFcn);
end


% Which model is used by default
if(~isfield(parameters,'model'))
    parameters.model = @model_Torord;
end


%% Running the model for the given number of beats.
CL = parameters.bcl;
time = cell(beats,1);
X = cell(beats, 1);
parameters.isFailed = 0; % We assume the simulation does not fail (which can change later).
for n=1:beats
    if (verbose)
        disp(['Beat = ' num2str(n)]);
    end
    
    if (maxTimePerBeat < Inf) % before each beat, restart the time counter, if time is to be measured
        timerVar = tic;
    end
    
    [time{n}, X{n}]=ode15s(parameters.model,[0 CL],X0,options,1,  cellType, ICaL_Multiplier, ...
        INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,...
        INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier,nao,cao,ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, apClamp, extraParams);
    
    
    if ~isequal(time{n}(end), parameters.bcl) % If simulation was killed prematurely in ode15sTimed, it is handled separately, setting dummy values.
        toc(timerVar);
        parameters.isFailed = 1; % It is noted the simulation was killed prematurely.
        try
            time(1:ignoreFirst) = [];
            X(1:ignoreFirst) = [];
            
        catch
            time = [];
            X = [];
        end
        
        return;
    end
    
    X0=X{n}(size(X{n},1),:);
end

%% Beats we're not interested in are deleted.
time(1:ignoreFirst) = [];
X(1:ignoreFirst) = [];


    %% A local function for limiting runtime - ignore if not using param.maxTimePerBeat.
    function [VALUE, ISTERMINAL, DIRECTION] = myEventFunction(T, Y, maxTimePerBeat, varargin)
    % varargin is there as our ode15s gets a lot of extra parameters, otherwise
    % the anonymous function handling this would complain and crash.

    %The event function stops when VALUE == 0 and
    %ISTERMINAL==1
    %a. Define the timeout in seconds
    TimeOut = maxTimePerBeat;
    %
    %b. The solver runs until this VALUE is negative (does not change the sign)
    VALUE = toc(timerVar)-TimeOut; % A nested function is used so we can access timerVar.
    %c. The function should terminate the execution, so
    ISTERMINAL = 1;
    %d. The direction does not matter
    DIRECTION = 0;
    end

end




