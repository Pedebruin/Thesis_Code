%{
    Since some of the simulation code such as the filter setup and the
    configuration is also necessary for the simulink files. The script is
    partitioned into several sub-scripts which can be run independantly:
        configure.m
        buildBeam.m
        buildFilters.m
%}
if exist('configuredExternally','var') == 1 && configuredExternally == true
    % Just leave the parameters that are in the workspace. 
else
    % Load all settings from the configuration file configure.m. All necessary
    % settings for this simulation, the model setup, filter setup and plotting
    % can be found there!

    clear
    close all
    set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);  
    
    % Configure the settings!
    [Beam, sBeam, modelSettings,...
            simulationSettings, plotSettings, simulinkSettings,...
            LO,KF,AKF,DKF,GDF,...
            Cd, configuredExternally] = configure(0);
end

now = char(datetime('now'));
fprintf([now,'-----------------------------\n'...
    'Setting up model...\n'])
startScript= tic;

stored = fileparts(which('MAIN.m'));
chdir(stored)
addpath('./Extras/')
addpath('./../c2000/')
addpath('./../c2000/Experiments/')

%% autoTune Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if simulationSettings.tuneFilters == true
    %{ 
        Find optimal values for filter parameters
        Use standard parameters as initial guess, then optimise further untill
        convergence. 
    %}
    options = optimset('Display','final',...
                       'MaxFunEvals',simulationSettings.maxFunEvals,...
                       'PlotFcn',@optimplotfval,...
                       'TolFun',eps);

    for i = 1:length(simulationSettings.filtersToTune)
        switch simulationSettings.filtersToTune
            case "AKF"
                disp(["Tuning AKF parameters... (Might take a while)"])
                filterTune = "AKF";
                AKF0 = [AKF.QTune,AKF.RTune];

                AKFopt = fminsearch(@(args)optfunc(args,filterTune),AKF0,options);

            %{
            Previous opt:
                QTuneAKF= 1.86e-11e-9;
                RTuneAKF= 2.26e4;
                QuTuned0AKF= 1e10;
            %}

            case "DKF"
                disp(["Tuning DKF parameters... (Might take a while)"])
                filterToTune= "DKF";
                DKF0 = [DKF.QTune,DKF.RTune];

                DKFopt = fminsearch(@(args)optfunc(args,filterToTune),DKF0,options);

                %{
                Previous opt:
                    QTuneDKF:3.0297e-13 
                    RTuneDKF: 1.1147e4
                    QuTuneDKF: 1e1
                %}

            case "GDF"
                disp(["Tuning GDF parameters... (Might take a while)"])
                filterToTune= "GDF";
                GDF0 = [GDF.QTune,GDF.RTune];

                GDFopt = fminsearch(@(args)optfunc(args,filterToTune),GDF0,options);

                %{
                Previous opt:
                    QTuneGDF: 1.90e-3
                    RTuneGDF: 1.65e14
                %}
        end
    end
  
        disp(["Bloody Tuned"])
else
    % Do nothing, and let the standard parameters take the wheel!
end

%% Build models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
For meta-analysis build multiple models with just different parameters. Put
them all in the modelMat array of models
%}

% Batch switch
if simulationSettings.batch == false
    simulationSettings.nLcurve = 1;
    simulationSettings.nDerivatives = 0;               % Reset to 0
    a = 1;
    b = 1;
    c = 1;
    d = 1;
else
    a = simulationSettings.nLcurve;             % First iteration dimension
    b = 1;                                      % Second iteration dimension
    c = simulationSettings.nDerivatives + 1;    % Third iteration dimension
    if simulationSettings.monteCarlo == true
        d = simulationSettings.iterations;          % MonteCarlo dimension
    else
        d = 1;
    end
end

% Initial error switch
if simulationSettings.initialError == false || strcmp(simulationSettings.data,"Real")
    simulationSettings.obsOffset = 0;
end

% MonteCarlo switch
if simulationSettings.monteCarlo == false
    simulationSettings.iterations = 1;
end

% Input measured noise covariance if real data
if strcmp(simulationSettings.data,"Real")
    load('noiseSample1.mat')

    laserNoise = squeeze(data{5}.Values.Data);
    laserNoise = laserNoise - mean(laserNoise);
    laserCov = cov(laserNoise');
    
    patchNoise = double(squeeze(data{2}.Values.Data));
    patchNoise = patchNoise -mean(patchNoise);
    patchCov = cov(patchNoise');

    accNoise = double(squeeze(data{1}.Values.Data))/simulationSettings.accSensitivity;
    accNoise = accNoise -mean(accNoise);
    accNoise = lowpass(accNoise,1e3,1e4);
    accCov = cov(accNoise');

    modelSettings.accCov = accCov;
    modelSettings.laserCov = laserCov;
    modelSettings.patchCov = patchCov;
end

% Set up modelMat (matrix with all models)
modelMat = model(simulationSettings.nLcurve,1,simulationSettings.nDerivatives+1);
simulationSettings.nModels = numel(modelMat);                                      

nPatches = length(modelSettings.patches);           % Number of patches
nsElementsP = modelSettings.nsElementsP;            % Number of elements per patch
nAcc = length(modelSettings.Acc);                   % Number of accelerometers

% Set batch vectors
LcurveQus = logspace(log10(simulationSettings.QuMin),log10(simulationSettings.QuMax),simulationSettings.nLcurve);

% Actually fill the modelMat matrix with appropriate models! (Varying batch
% settings)
l = 1; % Model counter

for i = 1:a % First iteration parameter
    for j = 1:b % Second iteration parameter
        for k = 1:c % Third iteration parameter
            %% Model setup
            mSettings = modelSettings;  % Copy of modelSettings to avoid overwriting the real thing. 
            mSettings.l = l;            % Also keep track of the iteration number for every model
            
            % True model for simulation
            SYS = model(); % create model object
            [SYS, mSettings] = buildBeam(mSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 
            
            % Save original modelSettings and beam
            originalBeam = copy(Beam);
            
            % Set modelling error
            if simulationSettings.modelError == true && strcmp(simulationSettings.data,"Simulated")
                Beam.rho = Beam.rho*mSettings.rhoError;                         % Beam mass error
                mSettings.mAcc = mSettings.mAcc*mSettings.mAccError;    % Accelerometer mass error
                mSettings.Acc = mSettings.Acc - mSettings.AccError;     % Acceleromter position error
            end
            
            % Erroneous model for the observers
            errSYS = model(); % create model object
            [errSYS,~] = buildBeam(mSettings,simulationSettings,plotSettings,Beam,sBeam,errSYS); % Fill model object
            
            if modelSettings.posMeasurement == true
                SYS.dsys_obs = errSYS.dsys_sim;    % Put only erronuous model in SYS and cut off laser measurement
            else
                SYS.dsys_obs = errSYS.dsys_sim(2:end,1:end);    % Put only erronuous model in SYS and cut off laser measurement
            end
            
            % Give summary of model
            if l == 1
                numEl = mSettings.numEl;
                Nmodes = mSettings.Nmodes;
                fprintf(['    # patches: %d \n'...
                         '    # accelerometers: %d \n'...
                         '    # elements: %d \n'...
                         '    # modes: %d \n'],nPatches,nAcc,numEl,Nmodes);
            end

            SYS.ny = size(SYS.dsys_sim,1);
            SYS.nu = size(SYS.dsys_sim,2);
            SYS.nq = size(SYS.dsys_sim.A,1);
            
            %% Filter setup
            % OverWrite derivativeOrder & QuTune
            AKF.nd = AKF.derivativeOrder;
            if simulationSettings.batch == true
                AKF.QuTune = LcurveQus(i);          % Change QuTune with i
                AKF.nd = 0;%k-1;                       % Change derivative order with k
            end

            % OverWrite derivativeOrder & QuTune
            DKF.nd = DKF.derivativeOrder;
            if simulationSettings.batch == true
                DKF.QuTune = LcurveQus(i);              % Change QuTune
                DKF.nd = k-1; % Change every loop!
            end

            SYS = buildFilters(SYS,LO,KF,AKF,DKF,GDF);
         
            %% CLosing remarks
            % Put SYS in modelMat matrix
            SYS.name = ['Model ', num2str(l)];
            SYS.number = l;
            modelMat(i,j,k) = SYS;

            % Reset modelSettings and beam
            Beam = originalBeam;
            l = l+1;
        end
    end
end

modelSettings = mSettings;

%% Set up all models to be simulated
if simulationSettings.simulate ==  true

    % Distubrance input generation (function at the bottom)
    switch simulationSettings.data
        case {"Real"}
            load(simulationSettings.dataset);
            t = data{1}.Values.Time;
            T = t(end);
            dt = mean(diff(t));
            
            Udist = -squeeze(data{4}.Values.Data)/simulinkSettings.matchingGain;
            
        case {"Simulated"}
            % Handy defenitions & Unpacking
            T = simulationSettings.T;
            dt = simulationSettings.dt;
            t = 0:dt:T;
            Udist = generateInput(simulationSettings);  % Disturbance input
            
            if simulationSettings.lowpassInput == true
                Udist = lowpass(Udist,simulationSettings.cutoffFrequency,1/simulationSettings.dt);
            end
            
            fprintf(['\n Setting up simulation...\n'...
            '    Noise: %d\n'...
            '    Offset: %1.1e\n'...
            '    Model Error: %d\n'],modelMat(1).simulationSettings.noise,modelMat(1).simulationSettings.obsOffset,modelMat(1).simulationSettings.modelError);
        
        otherwise 
            warning('SimulationSettings.data is not spelled correctly!')
    end

    fprintf(['    T: %.2f s \n'...
             '    dt: %.4f s\n'...
             '    Time steps: %d \n'],T, dt, length(t))

    if simulationSettings.batch == true
    fprintf(['    Batch:\n'...
             '        Number of models: %d \n'...
             '        Number of iterations per model: %d \n'...
             '        Number of total simulations: %d \n '],numel(modelMat),simulationSettings.iterations,simulationSettings.iterations*numel(modelMat));
    end
   
    %% Actual simulation loop!
    if simulationSettings.parallel == true
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool(simulationSettings.nWorkers);
        end
        if pool.NumWorkers == 1
            delete(pool)
            pool = parpool(simulationSettings.nWorkers);
        end
    else
        delete(gcp('nocreate'))     % Delete old one
        ps = parallel.Settings;
        ps.Pool.AutoCreate = false;
    end

    D = parallel.pool.DataQueue;

    % Simulate!!===========================================================
    l = 0;
    for i = 1:a
%     parfor i = 1:a                              
        for j = 1:b
            for k = 1:c
                SYS = modelMat(i,j,k); % Pick model to simulate
                for m = 1:d
                    warning('off','MATLAB:nearlySingularMatrix')
                    SYS = SYS.simulate(Udist,m); % Here the real magic happens!
                end
                SYS.simulationData.fit(:,1,SYS.simulationSettings.iterations+1) = mean(SYS.simulationData.fit(:,1,:),3);
                SYS.simulationData.fit(:,1,SYS.simulationSettings.iterations+2) = std(SYS.simulationData.fit(:,1,1:end-1),0,3);
                SYS.simulated = true;
                modelMat(i,j,k) = SYS; % put back in modelMat for safe keeping.   
             end
        end
    end
    
    % Save for later use if made with real dataset!
    if strcmp(modelMat(1).simulationSettings.data,"Real") && simulationSettings.saveResult == true
        save(['./dataSets/modelMat_',simulationSettings.dataset,],'modelMat')
    end
    
    %% Make a nice plot of the simulation!
    m = 1; % Only first iteration! (For illustrative purposes)
    if plotSettings.plot == true
        fprintf('\n Plotting simulation... \n')
        for i = 1:min(simulationSettings.nModels,20)
            simPlots = plotter(modelMat(i),m);
        end
    end
end
    
%% Evaluate batch info (Make L-curve plot)
if simulationSettings.nModels > 1
    lines = [];

    d = figure('Name','L-curve');
    hold on
    title('L-curve:');
    grid on
    xlabel('$$Q_u$$')
    ylabel('NRMSE [-]')
    patchAx = gca;
    set(patchAx,'xScale','log')
    set(patchAx,'yScale','log')
    ylim(gca,[1e-3,2])
    xlim(gca,[min(LcurveQus),max(LcurveQus)])

    for z = 1:length(modelMat(1,1).simulationSettings.observer)
        obs = simulationSettings.observer(z);
        if simulationSettings.monteCarlo == true
            means = zeros(a,b,c);
            sigmas = zeros(a,b,c);
            for i = 1:a
                for j = 1:b
                    for k = 1:c
                        means(i,j,k) = modelMat(i,j,k).simulationData.fit(z,1,end-1);
                        sigmas(i,j,k) = modelMat(i,j,k).simulationData.fit(z,1,end);
                    end
                end
            end
        else
            % RMSE
            fits = zeros(a,b,c);
            for i = 1:a
                for j = 1:b
                    for k = 1:c
                        fits(i,j,k) = modelMat(i,j,k).simulationData.fit(z,1,1);
                    end
                end
            end
        end
        
        switch obs 
            case {'AKF'}
                color = [0.4940 0.1840 0.5560];
            case {'DKF'}
                color = [0.3010 0.7450 0.9330];
        end
        
        if simulationSettings.monteCarlo == true
            meansk = means(:,1,k);
            sigmask = sigmas(:,1,k);
            switch obs
                case {'AKF'}
                    t = semilogx(patchAx,LcurveQus,meansk,'-o','Color',color);
                    patch(patchAx,[LcurveQus,fliplr(LcurveQus)],[meansk'+sigmask',fliplr(meansk'-sigmask')],color,'edgeAlpha',0);
                    lines = [lines,t];
                case {'DKF'}
                    t = semilogx(patchAx,LcurveQus,meansk,'-o','Color',color);
                    patch(patchAx,[LcurveQus,fliplr(LcurveQus)],[meansk'+sigmask',fliplr(meansk'-sigmask')],color,'edgeAlpha',0);
                    lines = [lines,t];
            end
            alpha(0.3)
        else
            t = loglog(patchAx,LcurveQus,fits(:,1,k),'-o','Color',color);
            lines = [lines,t];
        end
        movegui(d,"south")
    end
    legend(lines,simulationSettings.observer)
end

%% Also plot input and its derivatives
if plotSettings.inputPlot == true
    plotInput(modelMat(1),3);
end

% That's it, just some closing remarks. 
scriptElapsed = toc(startScript);
fprintf('DONE!!  in %2.2f Seconds! \n',scriptElapsed)

%%

%% FUNCTIONS!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input generation for impulse, step, harmonic and random signals
function [U] = generateInput(simulationSettings)


%{
Generates generic input sequences based on simulationSettings. 
Can generate: 
    Inpulse
    Step
    Harmonic
    Random
Adds these signals together in an imput signal U!
%}
    T = simulationSettings.T;
    dt = simulationSettings.dt;

    % Step input
    uStep = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.stepTime)                                % If step input is defined
        if simulationSettings.stepTime(1) == 0
            startSample = 1;
        else
            startSample = ceil(simulationSettings.stepTime(1)/dt);
        end

        uStep(startSample:end) = ones(1,T/dt-startSample+2)*simulationSettings.stepAmp;
        if length(simulationSettings.stepTime) > 1                          % If end time is defined
            stopSample = floor(simulationSettings.stepTime(2)/dt);
            uStep(stopSample:end) = zeros(1,T/dt-stopSample+2);
        end
    end

    % Inpulse input
    uImpulse = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.impulseTime)
        if simulationSettings.impulseTime == 0
            impulseSample = 1;
        else
            impulseSample = ceil(simulationSettings.impulseTime/dt);
        end
        uImpulse(impulseSample) = simulationSettings.impulseAmp;
    end

    % harmonic input
    uHarmonic = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.harmonicTime)
        if simulationSettings.harmonicTime == 0
            startSample = 1;
            t = 0:dt:T;
        else
            startSample = ceil(simulationSettings.harmonicTime(1)/dt);
            t = (simulationSettings.harmonicTime(1)-dt:dt:T)-simulationSettings.harmonicTime(1);
        end

        A = simulationSettings.harmonicAmp;
        f = simulationSettings.harmonicFreq;
        uHarmonic(startSample:end) = A*sin(2*pi*t*f);

        if length(simulationSettings.harmonicTime) > 1                          % If end time is defined
            stopSample = floor(simulationSettings.harmonicTime(2)/dt);
            uHarmonic(stopSample:end) = zeros(1,T/dt-stopSample+2);
        end
    end

    % Random Input
    uRand = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.randTime)
        if simulationSettings.harmonicTime == 0
            startSample = 1;
        else
            startSample = ceil(simulationSettings.randTime/dt);
        end
        a = simulationSettings.randInt(1);
        b = simulationSettings.randInt(2);
        uRand(startSample:end) = a + (b-a).*rand([1,T/dt-startSample+2]);
    end

    % Construct full signal!
    U = uStep + uImpulse + uHarmonic + uRand; 

    if size(0:dt:T) ~= length(U)
        error('Input length is not correct! Probably indexing issue :(')
    end
end

% Plot the input and its derivatives
function [] = plotInput(model,nd)

    simulationSettings = model.simulationSettings;
    
    if model.plotSettings.trueInputPlot == false
        % Generate custom input sequence to plot!
            simulationSettings.stepTime = [1,3];               % Location of input step ([start time], [endtime], [] )
                simulationSettings.stepAmp = 10;             % Step amplitude
            simulationSettings.impulseTime = [4];            % Location of input impulse ([time], [])
                simulationSettings.impulseAmp = 10;          % Inpulse amplitude
            simulationSettings.harmonicTime = [5,7];            % Harmonic input start time ([time], [])
                simulationSettings.harmonicFreq = 1;        % Frequency of sinusoidal input ([freq], [])
                simulationSettings.harmonicAmp = 10;         % Frequency input amplitude [Hz]
            simulationSettings.randTime = [8];               % random input start time ([time], [])
                simulationSettings.randInt = [-10,10];        % random input interval (uniformly distributed)
        
        simulationSettings.dt = 1e-3;
        simulationSettings.T = 10;
    end

    Udist = generateInput(simulationSettings);
    Udist = lowpass(Udist,1,1/simulationSettings.dt);

    t = 0:simulationSettings.dt:simulationSettings.T;

    U = zeros(nd+1,length(Udist));
    U(1,:) = Udist;

    for i = 1:nd+1
        U(i+1,1:end-1) = diff(U(i,:))/simulationSettings.dt;
    end
    
    derivatives = ["0^{th}" "1^{st}" "2^{nd}" "3^{rd}" "4^{th}" "5^{th}" "6^{th}"];
    units = ["$[N]$" "$[N/s]$" "$[N/s^2]$" "$[N/s^3]$"];

    figure('Name','Input & Derivatives')
    for i = 1:nd+1
        subplot(nd+1,1,i)
        hold on
        title(['$',char(derivatives(i)),'$', ' Derivative'])
        ylabel(units(i))
        plot(t,U(i,:))
    end
    xlabel('Time [s]')

end

% Optimise filters 
function [fit] = optfunc(args,filter)
%{
    args = [QTune, RTune, QuTuned0];
    Optimisation function that takes in set of parameters, simulates it to
    Udist and spits out the fit!
%}
    if any(args < 0)
        % Can't have negative covariances
        fit = NaN;
    else
        switch filter
            case "AKF"
                params = ["QTuneAKF","RTuneAKF"];%["QTuneAKF", "RTuneAKF","QuTuned0AKF"];
            case "DKF"
                params = ["QTuneDKF","RTuneDKF"];%["QTuneDKF", "RTuneDKF","QuTuned0DKF"];
            case "GDF"
                params = ["QTuneGDF", "RTuneGDF"];
        end
    
        [Beam, sBeam, modelSettings,...
        simulationSettings, plotSettings,simulinkSettings...
        ,LO,KF,AKF,DKF,GDF,~, ~] = configure(1,params,args,"Ssteel");
    
        AKF.nd = 0;
        DKF.nd = 0;
        modelSettings.l = 1;
        simulationSettings.observer = filter;
    
        % Input measured noise covariance if real data
        if strcmp(simulationSettings.data,"Real")
            load('noiseSample1.mat')
        
            laserNoise = squeeze(data{5}.Values.Data);
            laserNoise = laserNoise - mean(laserNoise);
            laserCov = cov(laserNoise');
            
            patchNoise = double(squeeze(data{2}.Values.Data));
            patchNoise = patchNoise -mean(patchNoise);
            patchCov = cov(patchNoise');
        
            accNoise = double(squeeze(data{1}.Values.Data))/simulationSettings.accSensitivity;
            accNoise = accNoise -mean(accNoise);
            accNoise = lowpass(accNoise,1e3,1e4);
            accCov = cov(accNoise');
        
            modelSettings.accCov = accCov;
            modelSettings.laserCov = laserCov;
            modelSettings.patchCov = patchCov;
        end
    
        SYS = model(); % create model object
        [SYS,~] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 
        
        if modelSettings.posMeasurement == true
            SYS.dsys_obs = SYS.dsys_sim;
        else
            SYS.dsys_obs = SYS.dsys_sim(2:end,1:end);    % Cut off laser measurement
        end
    
        SYS.ny = size(SYS.dsys_sim,1);
        SYS.nu = size(SYS.dsys_sim,2);
        SYS.nq = size(SYS.dsys_sim.A,1);
    
        SYS = buildFilters(SYS,LO,KF,AKF,DKF,GDF);    
    
        load(simulationSettings.dataset); % Load dataset to tune on
        Udist = -squeeze(data{4}.Values.Data)/simulinkSettings.matchingGain;
    
        warning('off','MATLAB:nearlySingularMatrix')
        SYS = SYS.simulate(Udist,1);
        warning('on','MATLAB:nearlySingularMatrix')
    
        if SYS.simulationData.broken == true
            fit = NaN;
        else
            fit = SYS.simulationData.fit;
            fprintf('P = [%10.4e, %10.4e]\n', args(1),args(2));
            disp(['Fit = ',num2str(fit(1))])
        end
    end
end