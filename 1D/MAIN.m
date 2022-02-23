clear
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

now = char(datetime('now'));
fprintf([now,'-----------------------------\n'...
    'Setting up model...\n'])
startScript= tic;

%% Parameters & Settings
% Basic settings
L = 150e-3;     % Beam length
b = 20e-3;      % Beam width
h = 0.1e-3;     % Beam thickness
patchL = 30e-3; % Patch length

% Model settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smart patches (Piezo)
    modelSettings.patches = [0];                          % location of start of patches (y/L)
    modelSettings.nsElementsP = 3;                          % Number of smart elements per patch
    modelSettings.LbElements = 0.1;                         % Preferred length of beam elements (y/L) (will change slightly)
    modelSettings.strainRate = false;                       % Measure strain rate instead of strain. 
    modelSettings.patchCov = -10;                    	    % True covariance of patch measurement (log10)
        modelSettings.maxPatchCov = 5;                      % Only when multiple models are simulated (log10)
        modelSettings.minPatchCov = -10;                    % Only when multiple models are simulated (log10)

% Accelerometers
    modelSettings.Acc = [1,0.75];                        % Location of accelerometers
    modelSettings.mAcc = 0.00127;                              % Mass of accelerometers
    modelSettings.accCov = -2;                             % True covariance of accelerometer measurement
        modelSettings.maxAccCov = 1;                      % Only when multiple models are simulated (log10)
        modelSettings.minAccCov = -10;                      % Only when multiple model sare simulated (log10)

% Strain gauges (TODO)

% Modelling 
    modelSettings.Nmodes = 5;                              % Number of modes to be modelled (Can't be larger then the amount of nodes)
    modelSettings.measurementHeight = 1;                    % Height of the measurement
    modelSettings.forceHeight = 1;                        % Height of the input force.
    
    modelSettings.wcov = 0;                              % Process noise covariance
    modelSettings.laserCov = 0;                          % Laser covariance                                          
    
    modelSettings.L = L;                                    
    modelSettings.b = b;
    
    modelSettings.rhoError = 1.05;                      % []% density error
    modelSettings.mAccError = 1.05;                     % []% accelerometer mass error
    modelSettings.AccError = 1e-3;                      % accelerometer position error
    
    modelSettings.c2dMethod = 'ZOH';
    
% simulationSettings (Settings for the simulation)%%%%%%%%%%%%%%%%%%%%%%%%%
simulationSettings.simulate = true;                     % Simulate at all?
    simulationSettings.waitBar = false;                      % Give waitbar and cancel option (VERY SLOW)
    
    % Turn on and off all errors!
    simulationSettings.noise = true;                   % Turn on or off all noise!
    simulationSettings.modelError = true;              % Turn on or off al model errors!
    simulationSettings.initialError = true;            % Turn on or off all initial state errors!
    
    % Time settings (Settings for the time and stepsize of the simulation)
    simulationSettings.dt = 1e-3;                       % Sampling time
    simulationSettings.T = 5;                           % Total simulation time
    
    % Input settings (Settings for the input that is used)
    simulationSettings.distInput = true;                   % Which input is the disturbance?
        simulationSettings.stepTime = [1,3];               % Location of input step ([start time], [endtime], [] )
            simulationSettings.stepAmp = 1e-2;             % Step amplitude
        simulationSettings.impulseTime = [];            % Location of input impulse ([time], [])
            simulationSettings.impulseAmp = 1e-2;          % Inpulse amplitude
        simulationSettings.harmonicTime = [];            % Harmonic input start time ([time], [])
            simulationSettings.harmonicFreq = 1;        % Frequency of sinusoidal input ([freq], [])
            simulationSettings.harmonicAmp = 1e-2;         % Frequency input amplitude [Hz]
        simulationSettings.randTime = [];               % random input start time ([time], [])
            simulationSettings.randInt = [-1,1]*1e-2;        % random input interval (uniformly distributed)
    simulationSettings.lowpassInput = false;
        simulationSettings.cutoffFrequency = 150;       % Hz

    %{
        Batch mode runs the script a lot of times to analyse meta
        results like the influence of different tuning parameters, or
        noise on performance. The script can vary 2 quantities, the
        'Patch' and 'Acc'elerometer covariance with batchMode (or
        both). It can also vary the amount of derivatives that ar
        emodelled. If nDerivatives is larger then 0, all lower
        derivatives are also modelled. Which increases simulation time
        significantly
    %}
    % Batch run settings
    simulationSettings.batch = false;                    % Run simulation in batch?? (A lot of times)
        simulationSettings.batchMode = 'None';            % 'Acc','Patch', 'Both' or 'None'
            simulationSettings.nPatchCov = 1 ;                  % How many data points between patchCov and minPatchCov?
            simulationSettings.nAccCov = 1;                    % How many data points between accCov and minaccCov?
        simulationSettings.nDerivatives = 0;               % Highest derivative order? (Also does all lower derivatives) 
      
    % Run each simulation multiple times for noise realisations?
    simulationSettings.monteCarlo = false;           % Run every simulation multiple times to get a mean and monteCarlo interval
        simulationSettings.iterations = 5;          % Amount of times every simulation is ran. 
   
    % Parallel computing settings
    simulationSettings.parallel = false;                % Run simulations in parallel?
        simulationSettings.nWorkers = 4;                % Number of cores to run this on?
    
    % Observer settings (Settings for the observers) 
    simulationSettings.observer = ["AKF"];    % ["MF" "LO" "KF" "AKF" "DKF" "GDF"] Does need to be in order
        simulationSettings.obsOffset = 1e-5;               % Initial state offset (we know its undeformed at the beginning, so probably 0); 

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
        AKF.QTune = 1e-6; % eye*QTune                  % Process noise covariance tuning
        AKF.RTune = 1; % R*RTune                        % Measurement noise covariance tuning
        AKF.P0 = 1e2;
        AKF.Pu0 = 1e5;

        AKF.QuTuned0 = 1e10;                             % For the zeroth derivative 
        AKF.QuTuned1 = 1e15;                             % For the first derivative
        AKF.QuTuned2 = 1e20;                             % For the second derivative

        %{
            Sinusoidal input A10f1: (Tuned)
                QuTuned0: 5e2
                QuTuned1: 1e4
                QuTuned2: 5e6; 

    	    Step input A10:  (Tuned)
                QuTuned0: 1e2
                QuTuned1: 1e6
                QuTuned2: 4e7
            
            Random input A10: (Not tuned)
                QuTuned0: 1e2
                QuTuned1: 0.5e2
                QuTuned2: 4e7
            
            Impulse input A10: (Not tuned)
                QuTuned0: 1e2
                QuTuned1: 0.5e2
                QuTuned2: 4e7
        %}

        % DKF settings
        DKF.derivativeOrder = 0;
        DKF.QTune = 1e-5; % eye*QTune
        DKF.RTune = 1; % R*RTune
        DKF.P0 = 1e-5;
        DKF.Pu0 = 1e-10;

        DKF.QuTuned0 = 1e-1;                                                                         
        DKF.QuTuned1 = 1e10;
        DKF.QuTuned2 = 1e20;

        %{
            Sinusoidal input A10f1: 
                QuTuned0: 
                QuTuned1: 
                QuTuned2:  

    	    Step input A10:
                QuTuned0: 1e1
                QuTuned1: 2.5e5
                QuTuned2: 
            
            Random input A10: 
                QuTuned0: 
                QuTuned1:
                QuTuned2:
            
            Impulse input A10:
                QuTuned0:
                QuTuned1:
                QuTuned2:
        %}

        % GDF settings
        GDF.QTune = 1e-6; % eye*QTune
        GDF.RTune = 1;  
        GDF.Pq0 = 1e-2;
        GDF.Pu0 = 1;
        GDF.Pqu0 = 0;

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

    plotSettings.inputPlot = false;                         % Plot the input and its derivatives?
        plotSettings.trueInputPlot = true;                  % Plot true or custom input sequence
                
    plotSettings.statePlot = true;                          % Plot the state evolutions
        plotSettings.states = 3;                            % First # states to be plotted

    plotSettings.modes = 0;                                 % Plot the mode shapes?? [number of modes]
        plotSettings.modeAmp = 5e-4;                        % Amplification factor for plot. 

    plotSettings.alphaMin = 0.3;                            % start of alpha interp

% beam element parameters (This is a beam element)
Beam = element('Beam');
Beam.h = h;                         % m thickness of beam
Beam.b = b;                         % m width of beam
Beam.E = 70E9;                      % E modulus of beam
Beam.mu = 0.334;                    % Poisson's ratio of beam
Beam.rho = 2710;                    % Mass density of beam
Beam.zeta = 0.01;                   % Modal damping coefficient of beam

% Smart beam parameters:
sBeam = copy(Beam); sBeam.name = 'sBeam';       
sBeam.L = patchL/modelSettings.nsElementsP;     % Length of smart beam element
sBeam.ph = 0.01e-3;                              % m 
sBeam.pb = b;                                   % m
sBeam.pE = 1e1;                      
sBeam.pmu = 0.334;
sBeam.prho = 1000;

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
    simulationSettings.batchMode = 'None';              % Choose correct covariances!
end

if simulationSettings.monteCarlo == false
    simulationSettings.iterations = 1;
end

% Set up modelMat (matrix with all models)
switch simulationSettings.batchMode
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
switch simulationSettings.batchMode
    case {'Patch'}
        patchCovariances = logspace(modelSettings.minPatchCov,modelSettings.maxPatchCov,simulationSettings.nPatchCov);
        patchCovariances = flip(patchCovariances);
    case {'Acc'}
        accCovariances = logspace(modelSettings.minAccCov,modelSettings.maxAccCov,simulationSettings.nAccCov);
        accCovariances = flip(accCovariances);
    case{'Both'}
        patchCovariances = logspace(modelSettings.minPatchCov,modelSettings.maxPatchCov,simulationSettings.nPatchCov);
        accCovariances = logspace(modelSettings.minAccCov,modelSettings.maxAccCov,simulationSettings.nAccCov);
    case{'None'}
        patchCovariances = 10^modelSettings.patchCov;
        accCovariances = 10^modelSettings.accCov;
end

% Initial error switch
if simulationSettings.initialError == false
    simulationSettings.obsOffset = 0;
end

% Actually fill the modelMat matrix with appropriate models! (Varying batch
% settings)
l = 1;
for i = 1:simulationSettings.nPatchCov
    for j = 1:simulationSettings.nAccCov
        for k = 1:simulationSettings.nDerivatives+1

        mSettings = modelSettings;  % Copy of modelSettings to avoid overwriting the real thing. 
       
        % Vary covariance of either patch or acc
        switch simulationSettings.batchMode
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
                mSettings.accCov = 10^modelSettings.accCov;
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


%% Set up all models to be simulated
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

    if simulationSettings.batch == true
    fprintf(['    Batch:\n'...
             '        Number of models: %d \n'...
             '        Number of iterations per model: %d \n'...
             '        Number of total simulations: %d \n '],numel(modelMat),simulationSettings.iterations,simulationSettings.iterations*numel(modelMat));
    end

    % Distubrance input generation (function at the bottom)
    Udist = generateInput(simulationSettings);  % Disturbance input
    
    if simulationSettings.lowpassInput == true
        Udist = lowpass(Udist,simulationSettings.cutoffFrequency,1/simulationSettings.dt);
    end

    % Set all models!
    l = 1;
    for i = 1:simulationSettings.nPatchCov
        for j = 1:simulationSettings.nAccCov
            for k = 1:simulationSettings.nDerivatives+1
                SYS = modelMat(i,j,k); % Pick model to set up
        
                SYS.ny = size(SYS.dsys_sim,1);
                SYS.nu = size(SYS.dsys_sim,2);
                SYS.nq = size(SYS.dsys_sim.A,1);
            
                % MF setup-------------------------------------------------------------
                if any(ismember(SYS.simulationSettings.observer,"MF"))
                    MF.Psi = pinv(SYS.dsys_obs.C);
                    SYS.MF = MF;
                else
                    MF = [];
                end
            
                % LO setup-------------------------------------------------------------
                if any(ismember(SYS.simulationSettings.observer,"LO"))
                    poles = eig(SYS.dsys_obs.A)*(1-LO.poleMovement);                   % Place poles in Discrete time (TODO)
                    LO.L = place(SYS.dsys_obs.A',SYS.dsys_obs.C',poles)'; 
            
                    SYS.LO = LO;
                    if l == 1
                    fprintf(['    Luenberger Observer-> \n'...
                             '        Pole speed: %3.1f%% increase \n'],LO.poleMovement*100)
                    end
                end
            
                % KF setup-------------------------------------------------------------
                if any(ismember(SYS.simulationSettings.observer,"KF"))
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
                    
                    SYS.KF = KF;
                    if l == 1
                        fprintf(['    Kalman Filter-> \n'...
                                 '        Stationary: %d \n'],KF.stationary)
                    end
                end
            
                % AKF setup------------------------------------------------------------
                AKF.nd = AKF.derivativeOrder;
                if any(ismember(SYS.simulationSettings.observer,"AKF"))
                    % q -> [q;qd;u] ->nq+nu*nd

                    % OverWrite derivativeOrder!
                    if simulationSettings.batch == true
                        AKF.nd = k-1; % Change every loop!
                    end

                    switch AKF.nd
                        case 0
                            % Use zeroth derivative tuning parameter
                            AKF.QuTune = AKF.QuTuned0;
                        case 1
                            % Use another tuning parameter
                            AKF.QuTune = AKF.QuTuned1;
                        case 2
                            % Use yet another tuning parameter
                            AKF.QuTune = AKF.QuTuned2;
                    end

                    AKF.S = [SYS.S(:,2:end); zeros(SYS.nu*(AKF.nd),SYS.ny-1)];

                    % Generate temporary ss model for input dynamics (higher order
                    % derivatives can then be modelled!
                        uA = [zeros(SYS.nu*AKF.nd,SYS.nu),eye(SYS.nu*AKF.nd);
                                zeros(SYS.nu,SYS.nu*(AKF.nd+1))];
                        uB = [zeros(AKF.nd*SYS.nu,SYS.nu);
                                eye(SYS.nu)];
                        uC = [eye(SYS.nu),zeros(SYS.nu,(AKF.nd)*SYS.nu)];
                        uD = [];
                        Uss = ss(uA,uB,uC,uD);
                        Uss = c2d(Uss,dt,SYS.modelSettings.c2dMethod); % Checked!
            
                    AKF.Dv = SYS.Dv(2:end,2:end);   % Snip off laser measurement
                    AKF.Bw = [SYS.Bw,zeros(SYS.nq,SYS.nu);
                        zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.B];
            
                    AKF.A = [SYS.dsys_obs.A, SYS.dsys_obs.B,zeros(SYS.nq,SYS.nu*AKF.nd);
                            zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.A];
                    AKF.C = [SYS.dsys_obs.C, SYS.dsys_obs.D, zeros(SYS.ny-1,SYS.nu*AKF.nd)];
                                              
                    AKF.R = SYS.R(2:end,2:end)*AKF.RTune;                                       % Do you trust your measurements?
                    AKF.Q = eye(size(SYS.Q,1))*AKF.QTune;                                           % Do you trust your model?
                    AKF.Qu = [1,zeros(1,nPatches);
                              zeros(nPatches,1),zeros(nPatches,nPatches)]*AKF.QuTune;                  % Input covariance!!
                    AKF.Qu(end-nPatches+1:end,end-nPatches+1:end) = eye(nPatches,nPatches)*1e-20;

                    AKF.Q = [AKF.Q,zeros(SYS.nq,SYS.nu); % Assemble augmented Q matrix
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

%                     AKFksys = minreal(ss(AKF.A,AKF.Bw,AKF.C,[],SYS.simulationSettings.dt));
%                     [~,~,P0_AKF] = kalman(AKFksys,AKF.Q,AKF.R,AKF.S);
%                     AKF.P0kalman = P0_AKF;

                    if l == 1
                        fprintf(['    Augmented Kalman Filter-> \n'...
                                 '        Stationary: %d \n'...
                                 '        QuTune: %1.1e \n'...
                                 '        nd: %d\n '],AKF.stationary,AKF.QuTune,AKF.nd)
                    end
                end
                SYS.AKF = AKF;
            
                % DKF setup------------------------------------------------------------
                DKF.nd = DKF.derivativeOrder;
                if any(ismember(simulationSettings.observer,"DKF"))
                    DKF.A = SYS.dsys_obs.A;
                    DKF.C = SYS.dsys_obs.C;
                    DKF.Bw = SYS.Bw;
                    

                    % OverWrite derivativeOrder!
                    if simulationSettings.batch == true
                        DKF.nd = k-1; % Change every loop!
                    end
                    
                    % Pick tuning parameter
                    switch DKF.nd
                        case 0
                            % Use zeroth derivative tuning parameter
                            DKF.QuTune = DKF.QuTuned0;
                        case 1
                            % Use another tuning parameter
                            DKF.QuTune = DKF.QuTuned1;
                        case 2
                            % Use yet another tuning parameter
                            DKF.QuTune = DKF.QuTuned2;
                    end
            
                        % Generate temporary ss model for input dynamics (higher order
                        % derivatives can then be modelled!
                        uA = [zeros(SYS.nu*DKF.nd,SYS.nu),eye(SYS.nu*DKF.nd);
                                zeros(SYS.nu,SYS.nu*(DKF.nd+1))];
                        uB = [zeros(DKF.nd*SYS.nu,SYS.nu);
                                eye(SYS.nu)];
                        uC = zeros(1,SYS.nu*(DKF.nd+1));
                        uC((0:SYS.nu-1)*(DKF.nd+1)+1) = ones(1,SYS.nu);
            
                        uD = [];
                        Uss = ss(uA,uB,uC,[]);
                        Uss = c2d(Uss,dt,modelSettings.c2dMethod);
                    
                    DKF.uA = Uss.A;
                    DKF.uB = Uss.B; 
                    DKF.uC = [eye(SYS.nu),zeros(SYS.nu,(DKF.nd)*SYS.nu)];
                    
                    DKF.D = SYS.dsys_obs.D;
                    DKF.B = SYS.dsys_obs.B;
                    DKF.Q = eye(size(SYS.Q,1))*DKF.QTune;
                    DKF.R = SYS.R(2:end,2:end)*DKF.RTune;                             % Snipp off laser measurement
                    DKF.Qu = eye(SYS.nu)*DKF.QuTune;
            
                    if l == 1
                        fprintf(['    Dual Kalman Filter-> \n'...
                         '        QuTune: %1.1e \n'],DKF.QuTune)
                    end
                end
                
                % Try and get the steady-state idare solution for both
                % systems, but only main system works. 
                % DKFkalmansys = ss(DKF.A,[DKF.B DKF.Bw],DKF.C,[DKF.D,zeros(SYS.ny-1,SYS.nq)],dt,'InputName','w','OutputName','y');
                % [~,~,~,~,P0_DKF] = kalman(DKFkalmansys,DKF.Q,DKF.R,[]);

                % DKFukalmansys = ss(DKF.uA,[DKF.uB,zeros(SYS.nu*(DKF.nd+1),SYS.nq)],DKF.D*DKF.uC,[zeros(SYS.ny-1,SYS.nu*(DKF.nd+1)), DKF.C],dt,'InputName','z','OutputName','y');
                % [~,~,~,~,Pu0_DKF] = kalman(DKFukalmansys,DKF.Qu,DKF.R,[]);

                SYS.DKF = DKF;
            
                % GDF setup------------------------------------------------------------
                if any(ismember(simulationSettings.observer,"GDF"))
                    GDF.A = SYS.dsys_obs.A;
                    GDF.B = SYS.dsys_obs.B;
                    GDF.C = SYS.dsys_obs.C;
                    GDF.D = SYS.dsys_obs.D;
            
                    GDF.Q = eye(size(SYS.Q,1))*GDF.QTune;
                    GDF.R = SYS.R(2:end,2:end)*GDF.RTune;                             % Snipp off laser measurement
            
                    
                    if l == 1
                        fprintf('    Giljins de Moor filter-> \n')
                    end
                end
                SYS.GDF = GDF;

                l = l+1;
                modelMat(i,j,k) = SYS; % put back in modelVec for safe keeping. 
            end
        end
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

    a = simulationSettings.nPatchCov;
    b = simulationSettings.nAccCov;
    c = simulationSettings.nDerivatives+1;
    d = SYS.simulationSettings.iterations;

    % Simulate!!===========================================================
    l = 0;
    for i = 1:a
    %parfor i = 1:a                              
        for j = 1:b
            for k = 1:c
                SYS = modelMat(i,j,k); % Pick model to simulate
                for m = 1:d 
                    SYS = SYS.simulate(Udist,m);
                end
                SYS.simulationData.fit(:,1,SYS.simulationSettings.iterations+1) = mean(SYS.simulationData.fit(:,1,:),3);
                SYS.simulationData.fit(:,1,SYS.simulationSettings.iterations+2) = std(SYS.simulationData.fit(:,1,1:end-1),0,3);
                SYS.simulated = 1;
                modelMat(i,j,k) = SYS; % put back in modelVec for safe keeping.   
             end
        end
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
    
%% Evaluate batch info
if simulationSettings.nModels > 1
    % Possibly save models (good practice after very long simulation)
    % save('patchAccCovBIG','modelMat')
    for z = 1:length(modelMat(1,1).simulationSettings.observer)
        if simulationSettings.monteCarlo == true
            means = zeros(simulationSettings.nPatchCov,simulationSettings.nAccCov,simulationSettings.nDerivatives);
            sigmas = zeros(simulationSettings.nPatchCov,simulationSettings.nAccCov,simulationSettings.nDerivatives);
            for i = 1:simulationSettings.nPatchCov
                for j = 1:simulationSettings.nAccCov
                    for k = 1:simulationSettings.nDerivatives+1
                        means(i,j,k) = modelMat(i,j,k).simulationData.fit(z,1,end-1);
                        sigmas(i,j,k) = modelMat(i,j,k).simulationData.fit(z,1,end);
                    end
                end
            end
        else
            % RMSE
            fits = zeros(simulationSettings.nPatchCov,simulationSettings.nAccCov,simulationSettings.nDerivatives);
            for i = 1:simulationSettings.nPatchCov
                for j = 1:simulationSettings.nAccCov
                    for k = 1:simulationSettings.nDerivatives+1
                        fits(i,j,k) = modelMat(i,j,k).simulationData.fit(z,1,1);
                    end
                end
            end
        end
        
        d = [];
        lines = [];
        colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
        for k = 1:simulationSettings.nDerivatives+1
        color = colors{k};
            switch simulationSettings.batchMode
                case {'Patch'}
                    if isempty(d)
                        d = figure('Name',join(['RMSE for increasing patch covariances',simulationSettings.observer(z)]));
                        hold on
                        title(join(['Patch covariance Vs estimation fit',simulationSettings.observer(z)]));
                        grid on
                        xlabel('Patch covariance $[V^2]$')
                        ylabel('NRMSE [-]')
                        patchAx = gca;
                        set(patchAx,'xScale','log')
                        set(patchAx,'yScale','log')
                        ylim(gca,[1e-3,1])
                        xlim(gca,[min(patchCovariances),max(patchCovariances)])
                    end

                    if simulationSettings.monteCarlo == true
                        meansk = means(:,1,k);
                        sigmask = sigmas(:,1,k);
                        t = semilogx(patchAx,patchCovariances,meansk,'-o','Color',color);
                        lines = [lines,t];
                        patch(patchAx,[patchCovariances,fliplr(patchCovariances)],[meansk'+sigmask',fliplr(meansk'-sigmask')],color,'edgeAlpha',0);
                        alpha(0.1)
                    else
                        loglog(patchAx,patchCovariances,fits(:,1,k),'-o','Color',color)
                    end
    
                case {'Acc'}
                    if isempty(d)
                        d = figure('Name','RMSE for increasing accelerometer covariances');
                        hold on
                        title 'Acc covariance Vs estimation fit'
                        grid on
                        xlabel('Accelerometer covariance $[(m/s^2)^2]$')
                        ylabel('NRMSE [-]')
                        accAx = gca;
                        set(accAx,'xScale','log')
                    end
                    semilogx(accAx,accCovariances,fits(1,:,k),'-o','Color',color)
    
                case {'Both'}
                    if isempty(d)
                    d = figure('Name','RMSE analysis');
                    hold on
                    title 'Acc covariance Vs Patch covariance'
                    grid on
                    xlabel('Patch covariance')
                    ylabel('Accelerometer covariance')
                    zlabel('RMSE of output')
                    surfAx = gca;
                    set(surfAx,'XScale','log')
                    set(surfAx,'YScale','log')
                    
                    end
                    surf(surfAx,patchCovariances,accCovariances,fits(:,:,k))
            end
            movegui(d,"south")
        end
        derivativeNames = ["0th order" "1st order" "2nd order" "3rd order" "4th order" "5th order" "6th order" "7th order" "8th order" "This is nuts"];
        legend(lines,derivativeNames(1:simulationSettings.nDerivatives+1))
    end
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
        uHarmonic(startSample:end) = A*sin(2*pi*t/f);

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