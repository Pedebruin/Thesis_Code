%{
This is a 1D Euler bernoulli simulation with piëzo patches

The beam model is generated based on Euler Bernoulli beam theory. The
piezoelec coupling is implemented as in the paper from 

Aktas, K. G.
Esen, I.
Doi:10.48084/etasr.3949

The author of this code is Pim de Bruin. 

The main model parameters are found at the top of the script.

The patches can be placed through modelSettings.patches. This is a vector
with the bottom location of each patch in natural coordinates. the length
of this vector determines the amount of patches. The size of the patches is 
determined through patchL. Same goes for the accelerometers and
modelSettings.Acc. 

The rest of the parameters are relatively self explanatory. If you have any
questions or find an error (very very probable) don't hesitate to send me a
message! @ P.E.deBruin@student.tudelft.nl

XXX Pim

This script requires:
    Symbolic math toolbox 
    Control system toolbox (For feedback command)
    Robust Control toolbox (For Hinf command)
    Statistics and machine learning toolbox (mvnrnd command)
%}

clear
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

now = char(datetime('now'));
fprintf([now,'-----------------------------\n'...
    'Setting up model...\n'])
startScript= tic;

%% Parameters & Settings
% Basic settings
L = 370e-3;     % Beam length
b = 40e-3;      % Beam width
h = 1e-3;       % Beam thickness
patchL = 50e-3; % Patch length

% Model settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smart patches (Piezo)
    modelSettings.patches = [0.8];                             % location of start of patches (y/L)
    modelSettings.nsElementsP = 3;                          % Number of smart elements per patch
    modelSettings.LbElements = 0.1;                         % Preferred length of beam elements (y/L) (will change slightly)
    modelSettings.strainRate = false;                        % Measure strain rate instead of strain. 
    modelSettings.patchCov = log10(1e2);                          % True covariance of patch measurement (log10)
        modelSettings.minPatchCov = -10;                   % Only when multiple models are simulated (log10)

% Accelerometers
    modelSettings.Acc = [0.95,0.8,0.2];                        % Location of accelerometers
    modelSettings.mAcc = 0.005;                              % Mass of accelerometers
    modelSettings.accCov = log10(1e2);                             % True covariance of accelerometer measurement
        modelSettings.minAccCov = -10;

% Strain gauges (TODO)

% Modelling 
    modelSettings.Nmodes = 5;                              % Number of modes to be modelled (Can't be larger then the amount of nodes)
    modelSettings.measurementHeight = 1;                    % Height of the measurement
    modelSettings.forceHeight = 0.1;                        % Height of the input force.
    
    modelSettings.wcov = 0;                              % Process noise covariance
    modelSettings.laserCov = 0;                          % Laser covariance                                          
    
    modelSettings.L = L;                                    
    modelSettings.b = b;

    modelSettings.rhoError = 1.05;                      % []% density error
    modelSettings.mAccError = 1.05;                     % []% accelerometer mass error
    modelSettings.AccError = 1e-4;                      % accelerometer position error

% simulationSettings (Settings for the simulation)%%%%%%%%%%%%%%%%%%%%%%%%%
simulationSettings.simulate = true;                     % Simulate at all?
    simulationSettings.metaAnalysis = 'Both';            % 'Acc','Patch', 'Both' or 'None'
    simulationSettings.waitBar = false;                      % Give waitbar and cancel option (VERY SLOW)
    
    % Turn on and off all errors!
    simulationSettings.noise = true;                   % Turn on or off all noise!
    simulationSettings.modelError = true;              % Turn on or off al model errors!
    simulationSettings.initialError = true;            % Turn on or off all initial state errors!
    
    % Time settings (Settings for the time and stepsize of the simulation)
    simulationSettings.dt = 1e-3;                       % Sampling time
    simulationSettings.T = 5;                           % Total simulation time
    
    % Input settings (Settings for the input that is used)
    simulationSettings.distInput = 1;                   % Which input is the disturbance?
        simulationSettings.stepTime = [];               % Location of input step ([start time], [endtime], [] )
            simulationSettings.stepAmp = 5;             % Step amplitude
        simulationSettings.impulseTime = [];            % Location of input impulse ([time], [])
            simulationSettings.impulseAmp = 10;          % Inpulse amplitude
        simulationSettings.harmonicTime = [0.1];            % Harmonic input start time ([time], [])
            simulationSettings.harmonicFreq = 1;        % Frequency of sinusoidal input ([freq], [])
            simulationSettings.harmonicAmp = 10;         % Frequency input amplitude [Hz]
        simulationSettings.randTime = [];               % random input start time ([time], [])
            simulationSettings.randInt = [-1,1];        % random input interval (uniformly distributed)
    
    % Batch run settings
    simulationSettings.batch = true;                    % Run simulation in batch?? (A lot of times)
        simulationSettings.nPatchCov = 10;                  % How many data points between patchCov and minPatchCov?
        simulationSettings.nAccCov = 10;                    % How many data points between accCov and minaccCov?
        simulationSettings.nDerivatives = 0;                % Number of derivatives?  
    
    % Observer settings (Settings for the observers) 
    simulationSettings.observer = ["AKF"];    % ["MF" "LO" "KF" "AKF" "DKF" "GDF"] Does need to be in order
        simulationSettings.obsOffset = 1e-5;                % Initial state offset (we know its undeformed at the beginning, so probably 0); 
        
        % MF settings

        % LO settings
        LO.poleMovement = 5e-3;                         % 0.5% faster poles

        % KF settings
        KF.stationary = false;	                        % Use stationary KF? 
        KF.QTune = 1;
        KF.RTune = 1;

        % AKF settings
        AKF.stationary = false;                         % Use stationary AKF? (DOES NOT WORK YET)
        AKF.derivativeOrder = 0;                        % Higher order derivative? (0:CP, 1:CV, 2:CA)
        AKF.QTune = 1;                                  % Process noise covariance tuning
        AKF.QuTune = 1e1;                               % Input sequence covariance tuning parameter
        AKF.RTune = 1;                                  % Measurement noise covariance tuning

        %{
            QuTune Patch:
            QuTune No Patch: 
        %}

        % DKF settings
        DKF.derivativeOrder = 0;                        % Does not work yet (TODO)
        DKF.QTune = 1;
        DKF.QuTune = 1;
        DKF.RTune = 1;

        % GDF settings
        GDF.QTune = 1;
        GDF.RTune = 1;

% plotSettings (Governs how the results are plotted!)%%%%%%%%%%%%%%%%%%%%%%
plotSettings.plot = true;                                   % Plot the simulation?
    plotSettings.plotNodes = true;                          % Plot the nodes
        plotSettings.nodeNumbers = true;                    % Plot the node numbers
    plotSettings.elementNumbers = false;                    % Plot the element numbers
    plotSettings.piezos = true;                             % Plot the piezo elements
        plotSettings.piezoNumbers = true;                   % Give them numbers
    plotSettings.sensor = true;                             % Plot the sensor in beam plot (red line)
    plotSettings.inputForce = true;                         % Plot force location in beam plot
    plotSettings.inputSequence = true;                      % Plot input force sequence?
    plotSettings.accelerometers = true;                     % Plot accelerometers?
        plotSettings.accNumbers = true;                     % Plot acceleromter numbers?

    plotSettings.statePlot = true;                          % Plot the state evolutions
        plotSettings.states = 5;                            % First # states to be plotted

    plotSettings.modes = 0;                                 % Plot the mode shapes?? [number of modes]
        plotSettings.modeAmp = 5e-3;                        % Amplification factor for plot. 

    plotSettings.alphaMin = 0.3;                            % start of alpha interp

% beam element parameters (This is a beam element)
Beam = element('Beam');
Beam.h = h;                         % m thickness of beam
Beam.b = b;                         % m width of beam
Beam.E = 70E9;                      % E modulus of beam
Beam.mu = 0.334;                    % Poisson's ratio of beam
Beam.rho = 2710;                    % Mass density of beam
Beam.zeta = 0.01;                   % Modal damping coefficient of beam

% Smart beam parameters :PI P-876 DuraAct (This is a 'smart' element)
sBeam = copy(Beam); sBeam.name = 'sBeam';       
sBeam.L = patchL/modelSettings.nsElementsP;     % Length of smart beam element
sBeam.ph = 0.5e-3;                              % m 
sBeam.pb = b;                                   % m
sBeam.pE = 5.2e10;                      
sBeam.pmu = 0.334;
sBeam.prho = 7800;

% strain element parameters (TODO)

%% Build beam model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
For meta-analysis build multiple models with just different parameters. Put
them all in the modelVec vector of models
%}

if simulationSettings.batch == false
    % Just to be sure, built a master switch to keep from accidentally
    % creating too much simulations. 
    simulationSettings.nPatchCov = 1;                  % Reset to 1
    simulationSettings.nAccCov = 1;                    % Reset to 1
    simulationSettings.nDerivatives = 0;               % Reset to 0
end

switch simulationSettings.metaAnalysis
    case {'Patch'}
        modelMat = model(simulationSettings.nPatchCov,1,simulationSettings.nDerivatives+1);
        simulationSettings.nAccCov = 1;
    case {'Acc'}
        modelMat = model(1,simulationSettings.nAccCov,simulationSettings.nDerivatives+1);
        simulationSettings.nPatchCov = 1;
    case {'Both'}
        modelMat = model(simulationSettings.nPatchCov,simulationSettings.nAccCov,simulationSettings.nDerivatives+1);
    case {'None'}
        modelMat = model(1,1,simulationSettings.nDerivatives+1);
        simulationSettings.nAccCov = 1;
        simulationSettings.nPatchCov = 1;
end

simulationSettings.nModels = numel(modelMat);                                      

nPatches = length(modelSettings.patches);           % Number of patches
nsElementsP = modelSettings.nsElementsP;            % Number of elements per patch
nAcc = length(modelSettings.Acc);                   % Number of accelerometers

% generate meta analysis vector
switch simulationSettings.metaAnalysis
    case {'Patch'}
        patchCovariances = logspace(modelSettings.minPatchCov,modelSettings.patchCov,simulationSettings.nPatchCov);
        patchCovariances = flip(patchCovariances);
    case {'Acc'}
        accCovariances = logspace(modelSettings.minAccCov,modelSettings.accCov,simulationSettings.nAccCov);
        accCovariances = flip(accCovariances);
    case{'Both'}
        patchCovariances = logspace(modelSettings.minPatchCov,modelSettings.patchCov,simulationSettings.nPatchCov);
        accCovariances = logspace(modelSettings.minAccCov,modelSettings.accCov,simulationSettings.nAccCov);
    case{'None'}
        patchCovariances = 10^modelSettings.patchCov;
        accCovariances = 10^modelSettings.accCov;
end

% Initial error switch
if simulationSettings.initialError == false
    simulationSettings.obsOffset = 0;
end

l = 1;
for i = 1:simulationSettings.nPatchCov
    for j = 1:simulationSettings.nAccCov
        for k = 1:simulationSettings.nDerivatives+1

        mSettings = modelSettings;  % Copy of modelSettings to avoid overwriting the real thing. 
       
        % Vary covariance of either patch or acc
        switch simulationSettings.metaAnalysis
            case {'Patch'}
                mSettings.patchCov = patchCovariances(i);
                mSettings.accCov = 10^mSettings.accCov;
            case {'Acc'}
                mSettings.accCov = accCovariances(j);
                mSettings.patchCov = 10^mSettings.patchCov;
            case {'Both'}
                mSettings.patchCov = patchCovariances(i);
                mSettings.accCov = accCovariances(j);
            case {'None'}
                mSettings.patchCov = 10^modelSettings.patchCov;
                mSettings.accCov = 10^modelSettings.patchCov;
        end
    
        % True model for simulation
        SYS = model(); % create model object
        [SYS, mSettings] = buildBeam(mSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 
    
        % Save original modelSettings and beam
        originalBeam = copy(Beam);
    
        % Set modelling error
        if simulationSettings.modelError == true
            Beam.rho = Beam.rho*mSettings.rhoError;                         % Beam mass error
            mSettings.mAcc = mSettings.mAcc*mSettings.mAccError;    % Accelerometer mass error
            mSettings.Acc = mSettings.Acc - mSettings.AccError;     % Acceleromter position error
        end
    
        errSYS = model(); % create model object
        [errSYS,~] = buildBeam(mSettings,simulationSettings,plotSettings,Beam,sBeam,errSYS); % Fill model object
        SYS.dsys_obs = errSYS.dsys_sim(2:end,1:end);    % Put only erronuous model in SYS and cut off laser measurement
    
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


% Give summary of model
numEl = modelSettings.numEl;
Nmodes = modelSettings.Nmodes;
fprintf(['    # patches: %d \n'...
         '    # accelerometers: %d \n'...
         '    # elements: %d \n'...
         '    # modes: %d \n'],nPatches,nAcc,numEl,Nmodes);


%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulates the system and makes nice plots and stuff. It is not
% really necessary and not the point of this script. But worked wonders when
% debugging the code. 

if simulationSettings.simulate ==  true
    fprintf(['\n Setting up simulation...\n'...
            '    Noise: %d\n'...
            '    Offset: %1.1e\n'...
            '    Model Error: %d\n'],modelMat(1).simulationSettings.noise,modelMat(1).simulationSettings.obsOffset,modelMat(1).simulationSettings.modelError);

    % Handy defenitions & Unpacking
    T = simulationSettings.T;
    dt = simulationSettings.dt;
    t = 0:dt:T;

    fprintf(['    T: %.2f s \n'...
    '    dt: %.4f s\n'...
    '    Time steps: %d \n'],simulationSettings.T, simulationSettings.dt, length(t))

    % Distubrance input generation (function at the bottom)
    Udist = generateInput(simulationSettings);  % Disturbance input

    % Set up waitbar if necessary
    if simulationSettings.waitBar == true
        f = waitbar(0,'1','Name','Running simulation loop...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(f,'canceling',0);
    end
    
    % Set up filters and simulate multiple times!
    l = 1;
    for i = 1:simulationSettings.nPatchCov
        for j = 1:simulationSettings.nPatchCov
            for k = 1:simulationSettings.nDerivatives+1
                SYS = modelMat(i,j,k); % Pick model to simulate
        
                SYS.ny = size(SYS.dsys_sim,1);
                SYS.nu = size(SYS.dsys_sim,2);
                SYS.nq = size(SYS.dsys_sim.A,1);
            
                % MF setup-------------------------------------------------------------
                if any(ismember(simulationSettings.observer,"MF"))
                    MF.Psi = pinv(SYS.dsys_obs.C);
                else
                    MF = [];
                end
            
                % LO setup-------------------------------------------------------------
                if any(ismember(simulationSettings.observer,"LO"))
                    poles = eig(SYS.dsys_obs.A)*(1-LO.poleMovement);                   % Place poles in Discrete time (TODO)
                    LO.L = place(SYS.dsys_obs.A',SYS.dsys_obs.C',poles)'; 
            
                    if l == 1
                    fprintf(['    Luenberger Observer-> \n'...
                             '        Pole speed: %3.1f%% increase \n'],LO.poleMovement*100)
                    end
                end
            
                % KF setup-------------------------------------------------------------
                if any(ismember(simulationSettings.observer,"KF"))
                    KF.Q = SYS.Q*KF.QTune;
                    KF.R = SYS.R(2:end,2:end)*KF.RTune;                          % Snip off laser measurement
                    KF.S = SYS.S(:,2:end);
                    KF.Bw = SYS.Bw;
                    KF.Dv = SYS.Dv(2:end,2:end); 
                      
                    if KF.stationary == true
                        % Find stationary kalman gain by solving discrete algebraic ricatti equation
                        % (p.162 M.Verhaegen
                        [~,Kt,KF.poles] = idare(SYS.dsys_obs.A, SYS.dsys_obs.C', KF.Q, KF.R, KF.S, []);
                        KF.K = Kt';
                    else
                        % No additional setup for non-stationary kalman filter.
                    end
                    
                    if l == 1
                        fprintf(['    Kalman Filter-> \n'...
                                 '        Stationary: %d \n'],KF.stationary)
                    end
                end
            
                % AKF setup------------------------------------------------------------
                AKF.nd = AKF.derivativeOrder;
                if any(ismember(simulationSettings.observer,"AKF"))
                    % q -> [q;qd;u] ->nq+nu*nd
                    AKF.nd = k-1; % Change every loop!

                    AKF.S = [SYS.S(:,2:end); zeros(SYS.nu*(AKF.nd+1),SYS.ny-1)];
            
                    % Generate temporary ss model for input dynamics (higher order
                    % derivatives can then be modelled!
                        uA = [zeros(SYS.nu*AKF.nd,SYS.nu),eye(SYS.nu*AKF.nd);
                                zeros(SYS.nu,SYS.nu*(AKF.nd+1))];
                        uB = [zeros(AKF.nd*SYS.nu,SYS.nu);
                                eye(SYS.nu)];
                        uC = zeros(1,SYS.nu*AKF.nd);
                        uD = [];
                        Uss = ss(uA,uB,[],[]);
                        Uss = c2d(Uss,dt,modelSettings.c2dMethod);
            
                    AKF.Dv = SYS.Dv(2:end,2:end);
                    AKF.Bw = [SYS.Bw,zeros(SYS.nq,SYS.nu);
                        zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.B];
            
                    AKF.A = [SYS.dsys_obs.A, SYS.dsys_obs.B,zeros(SYS.nq,SYS.nu*AKF.nd);
                            zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.A];
                    AKF.C = [SYS.dsys_obs.C, SYS.dsys_obs.D, zeros(SYS.ny-1,SYS.nu*AKF.nd)];
                                              
                    AKF.R = SYS.R(2:end,2:end)*AKF.RTune;                                       % Do you trust your measurements?
                    AKF.Q = SYS.Q*AKF.QTune;                                           % Do you trust your model?
                    AKF.Qu = [1,zeros(1,nPatches);
                              zeros(nPatches,1),zeros(nPatches,nPatches)]*AKF.QuTune;                  % Input covariance!!
            
                    AKF.Q = [SYS.Q,zeros(SYS.nq,SYS.nu); % Assemble augmented Q matrix
                        zeros(SYS.nu,SYS.nq),AKF.Qu];
                
                    if AKF.stationary == true
                        % Find stationary kalman gain by solving discrete algebraic ricatti equation
                        % (p.162 M.Verhaegen)
                        if AKF.nd > 0 || AKF.nd > 1
                            error('Ricatti does not work for this augmented system.. (Already in backlog, working on it)')
                        end
                        [~,Kt,AKF.poles,INFO] = idare(AKF.A, AKF.C', AKF.Q, AKF.R, AKF.S, []);
                        AKF.K = Kt';
                    else
                        % No additional setup for non-stationary augmented kalman filter.
                    end
            
                    if l == 1
                        fprintf(['    Augmented Kalman Filter-> \n'...
                                 '        Stationary: %d \n'...
                                 '        QuTune: %1.1e \n'],AKF.stationary,AKF.QuTune)
                    end
                end
            
                % DKF setup------------------------------------------------------------
                if any(ismember(simulationSettings.observer,"DKF"))
                    nd_DKF = DKF.derivativeOrder;
                    DKF.A = SYS.dsys_obs.A;
                    DKF.B = SYS.dsys_obs.B;
                    DKF.C = SYS.dsys_obs.C;
                    DKF.D = SYS.dsys_obs.D;
            
                        % Generate temporary ss model for input dynamics (higher order
                        % derivatives can then be modelled!
                        uA = [zeros(SYS.nu*nd_DKF,SYS.nu),eye(SYS.nu*nd_DKF);
                                zeros(SYS.nu,SYS.nu*(nd_DKF+1))];
                        uB = [zeros(nd_DKF*SYS.nu,SYS.nu);
                                eye(SYS.nu)];
                        uC = zeros(1,SYS.nu*(nd_DKF+1));
                        uC((0:SYS.nu-1)*(nd_DKF+1)+1) = ones(1,SYS.nu);
            
                        uD = [];
                        Uss = ss(uA,uB,uC,[]);
                        Uss = c2d(Uss,dt,modelSettings.c2dMethod);
                    
                    DKF.uA = Uss.A;
                    DKF.uB = Uss.B; 
                    DKF.uC = Uss.C; % IMPLEMENT THIS IN SIMULATION LOOP! (So does not work yet)
            
                    DKF.Q = SYS.Q*DKF.QTune;
                    DKF.R = SYS.R(2:end,2:end)*DKF.RTune;                             % Snipp off laser measurement
                    DKF.Qu = eye(SYS.nu)*DKF.QuTune;
            
                    if l == 1
                        fprintf(['    Dual Kalman Filter-> \n'...
                         '        QuTune: %1.1e \n'],DKF.QuTune)
                    end
                end
            
                % GDF setup------------------------------------------------------------
                if any(ismember(simulationSettings.observer,"GDF"))
                    GDF.A = SYS.dsys_obs.A;
                    GDF.B = SYS.dsys_obs.B;
                    GDF.C = SYS.dsys_obs.C;
                    GDF.D = SYS.dsys_obs.D;
            
                    GDF.Q = SYS.Q*GDF.QTune;
                    GDF.R = SYS.R(2:end,2:end)*GDF.RTune;                             % Snipp off laser measurement
            
                    if l == 1
                        fprintf('    Giljins de Moor filter-> \n')
                    end
                end
        
                % Simulate!!===========================================================
                SYS = SYS.simulate(MF,LO,KF,AKF,DKF,GDF,Udist);
                % =====================================================================
                
                l = l+1;
                modelMat(i,j,k) = SYS; % put back in modelVec for safe keeping. 
            end
        end
    end
    
    % Make a nice plot of the simulation!
    if plotSettings.plot == true
        fprintf('\n Plotting simulation... \n')
        for i = 1:min(simulationSettings.nModels,15)
            simPlots = plotter(modelMat(i));
        end
    end
end
    

%% Evaluate batch info
if simulationSettings.nModels > 1
    % Possibly save models
    % save('patchAccCovBIG','modelMat')

    % RMSE
    fits = zeros(simulationSettings.nPatchCov,simulationSettings.nAccCov);
    y_true =  modelMat(1).simulationData.yfull(1,:); % They all have the same true simulation

    for i = 1:simulationSettings.nPatchCov
        for j = 1:simulationSettings.nAccCov
            y_AKF = modelMat(i,j).simulationData.yfull_AKF(1,:);
            fits(i,j) = goodnessOfFit(y_AKF',y_true','NRMSE');
        end
    end
    
    switch simulationSettings.metaAnalysis
        case {'Patch'}
            d = figure('Name','RMSE for increasing patch covariances');
            semilogx(patchCovariances,fits,'-o')
            hold on
            title 'Patch covariance Vs estimation fit'
            grid on
            xlabel('Patch covariance $[V^2]$')
            ylabel('NRMSE [-]')
        case {'Acc'}
            d = figure('Name','RMSE for increasing accelerometer covariances');
            semilogx(accCovariances,fits,'-o')
            hold on
            title 'Acc covariance Vs estimation fit'
            grid on
            xlabel('Accelerometer covariance $[(m/s^2)^2]$')
            ylabel('NRMSE [-]')
        case {'Both'}
            d = figure('Name','RMSE analysis');
            surf(patchCovariances,accCovariances,fits)
            hold on
            title 'Acc covariance Vs Patch covariance'
            grid on
            xlabel('Patch covariance')
            ylabel('Accelerometer covariance')
            zlabel('RMSE of output')

            set(gca,'XScale','log')
            set(gca,'YScale','log')
    end
    movegui(d,"south")
end

scriptElapsed = toc(startScript);

fprintf('DONE!!  in %2.2f Seconds! \n',scriptElapsed)

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
            startSample = ceil(simulationSettings.harmonicTime/dt);
            t = (simulationSettings.harmonicTime-dt:dt:T)-simulationSettings.harmonicTime;
        end
        
        A = simulationSettings.harmonicAmp;
        f = simulationSettings.harmonicFreq;
        uHarmonic(startSample:end) = A*sin(2*pi*t/f);
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