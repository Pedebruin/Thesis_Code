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
h = 1e-3;     % Beam thickness

% Model settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beam element parameters (This is a beam element)
    Beam = element('Beam');
    Beam.h = h;                         % m thickness of beam
    Beam.b = b;                         % m width of beam
    Beam.E = 70E9;                      % E modulus of beam
    Beam.mu = 0.334;                    % Poisson's ratio of beam
    Beam.rho = 2710;                    % Mass density of beam
    Beam.zeta = 0.01;                   % Modal damping coefficient of beam 

% Smart patches (Piezo or strain)
    modelSettings.patches = [0];                                % location of start of patches (y/L)
    modelSettings.gauges = "resistive";                         % "resistive" or "piezo" gauges? 
    modelSettings.nsElementsP = 3;                              % Number of smart elements per patch
    modelSettings.strainRate = false;                           % Measure strain rate instead of strain. (For current amplifiers)
    modelSettings.posMeasurement = false;                        % Use position measurement for the observers?

    sBeam = copy(Beam); sBeam.name = 'sBeam';   
    switch modelSettings.gauges
        case {'piezo'}                                          % When using piezo patches
            patchL = 30e-3;                                     % Patch length
            sBeam.ph = 0.5e-3;                                  % m 
            sBeam.pb = b;                                       % m
            sBeam.pE = 5.2e10;                      
            sBeam.pmu = 0.334;
            sBeam.prho = 7800;
            modelSettings.patchCov = 1e-5;                    	% True covariance of patch measurement (log10)
        case {'resistive'}                                      % When using resistive strain gauges
            patchL = 30e-3;                                     % Patch length
            sBeam.ph = 0.1e-3;                                  % m 
            sBeam.pb = b;                                       % m
            sBeam.pE = 10;                      
            sBeam.pmu = 0.334;
            sBeam.prho = 1000;
            modelSettings.patchCov = 1e-13;                       % True covariance of patch measurement (log10)
    end
    sBeam.L = patchL/modelSettings.nsElementsP;                 % Length of smart beam element

% Accelerometers
    modelSettings.Acc = [1,0.75];                        	    % Location of accelerometers
    modelSettings.mAcc = 0.00127;                               % Mass of accelerometers
    modelSettings.accCov = 1e-2;                                  % True covariance of accelerometer measurement

% General Modelling 
    modelSettings.LbElements = 0.1;                             % Preferred length of beam elements (y/L) (will change slightly)
    modelSettings.Nmodes = 5;                                   % Number of modes to be modelled (Can't be larger then the amount of nodes)
    modelSettings.measurementHeight = 1;                        % Height of the measurement
    modelSettings.forceHeight = 1;                              % Height of the input force.
    
    modelSettings.wcov = 0;                              % Process noise covariance
    modelSettings.laserCov = 1e-12;                          % Laser covariance                                          
    
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
            simulationSettings.stepAmp = 1;             % Step amplitude
        simulationSettings.impulseTime = [];            % Location of input impulse ([time], [])
            simulationSettings.impulseAmp = 1;          % Inpulse amplitude
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
    simulationSettings.batch = false;                   % Run simulation in batch?? (A lot of times)
        simulationSettings.nDerivatives = 0;            % Highest derivative order? (Also does all lower derivatives) 
        simulationSettings.QuMin = 1e1;                 % Loop from
        simulationSettings.QuMax = 1e12;                % Loop to
        simulationSettings.nLcurve = 25;                 % In .. steps

    % Run each simulation multiple times for noise realisations?
    simulationSettings.monteCarlo = false;           % Run every simulation multiple times to get a mean and monteCarlo interval
        simulationSettings.iterations = 5;          % Amount of times every simulation is ran. 
   
    % Parallel computing settings
    simulationSettings.parallel = false;                % Run simulations in parallel?
        simulationSettings.nWorkers = 4;                % Number of cores to run this on?
    
    % Observer settings (Settings for the observers) %%%%%%%%%%%%%%%%%%%%%%
    simulationSettings.observer = ["AKF" "DKF"];    % ["MF" "LO" "KF" "AKF" "DKF" "GDF"] Does need to be in order
        simulationSettings.obsOffset = 1e-10;               % Initial state offset (we know its undeformed at the beginning, so probably 0); 
    
        % LO settings
        LO.poleMovement = 5e-3;                         % 0.5% faster poles

        % KF settings
        KF.stationary = false;	                        % Use stationary KF? 
        KF.QTune = 1e-10;
        KF.RTune = 1e10;

        % AKF settings
        AKF.stationary = false;                         % Use stationary AKF? (DOES NOT WORK YET)
        AKF.derivativeOrder = 0;                        % % Derivative order when not in batch mode (0:CP, 1:CV, 2:CA)
        AKF.QTune = 1e-10; % eye*QTune                  % Process noise covariance tuning
        AKF.RTune = 1e10; % R*RTune                        % Measurement noise covariance tuning
        AKF.P0 = 0;
        AKF.Pu0 = 0;

        AKF.QuTuned0 = 5e10;                             % For the zeroth derivative 
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
        DKF.derivativeOrder = 0;                                                                            % Derivative order when not in batch mode
        DKF.QTune = 1e-30; % eye*QTune
        DKF.RTune = 1e10; % R*RTune
        DKF.P0 = 0;
        DKF.Pu0 = 0;

        DKF.QuTuned0 = 5e10;                                                                         
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
        GDF.QTune = 1e-1; % eye*QTune
        GDF.RTune = 1;  
        GDF.Pq0 = 0;
        GDF.Pu0 = 0;
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

%% Build models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
For meta-analysis build multiple models with just different parameters. Put
them all in the modelMat arra of models
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
if simulationSettings.initialError == false
    simulationSettings.obsOffset = 0;
end

% MonteCarlo switch
if simulationSettings.monteCarlo == false
    simulationSettings.iterations = 1;
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
                if modelSettings.posMeasurement == true
                    KF.R = SYS.R*KF.RTune;
                    KF.S = SYS.S;
                    KF.Dv = SYS.Dv;
                else
                    KF.R = SYS.R(2:end,2:end)*KF.RTune;                          % Snip off laser measurement
                    KF.S = SYS.S(:,2:end);
                    KF.Dv = SYS.Dv(2:end,2:end); 
                end
                KF.Bw = SYS.Bw;
                
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
    
                % OverWrite derivativeOrder % QuTune
                if simulationSettings.batch == true
                    AKF.QuTune = LcurveQus(i);          % Change QuTune with i
                    AKF.nd = k-1;                       % Change derivative order with k
                end

    
                % Generate temporary ss model for input dynamics (higher order
                % derivatives can then be modelled!
                    uA = [zeros(SYS.nu*AKF.nd,SYS.nu),eye(SYS.nu*AKF.nd);
                            zeros(SYS.nu,SYS.nu*(AKF.nd+1))];
                    uB = [zeros(AKF.nd*SYS.nu,SYS.nu);
                            eye(SYS.nu)];
                    uC = [eye(SYS.nu),zeros(SYS.nu,(AKF.nd)*SYS.nu)];
                    uD = [];
                    Uss = ss(uA,uB,uC,uD);
                    Uss = c2d(Uss,simulationSettings.dt,mSettings.c2dMethod); % Checked!
        
                if modelSettings.posMeasurement == true % Depending on settings, snip off laser measurement. 
                    AKF.R = SYS.R*AKF.RTune; 
                    AKF.S = [SYS.S; zeros(SYS.nu*(AKF.nd),SYS.ny-1)];
                    AKF.Dv = SYS.Dv;
                else
                    AKF.R = SYS.R(2:end,2:end)*AKF.RTune;                                       % Do you trust your measurements?
                    AKF.S = [SYS.S(:,2:end); zeros(SYS.nu*(AKF.nd),SYS.ny-1)];
                    AKF.Dv = SYS.Dv(2:end,2:end);
                end

                AKF.Bw = [SYS.Bw,zeros(SYS.nq,SYS.nu);
                    zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.B];
                AKF.A = [SYS.dsys_obs.A, SYS.dsys_obs.B,zeros(SYS.nq,SYS.nu*AKF.nd);
                    zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.A];
                AKF.C = [SYS.dsys_obs.C, SYS.dsys_obs.D, zeros(SYS.ny-1,SYS.nu*AKF.nd)];
                                          
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
    
                if l == 1
                    fprintf(['    Augmented Kalman Filter-> \n'...
                             '        Stationary: %d \n'...
                             '        nd: %d\n '],AKF.stationary,AKF.nd)
                end
            end
            SYS.AKF = AKF;
        
            % DKF setup------------------------------------------------------------
            DKF.nd = DKF.derivativeOrder;
            if any(ismember(simulationSettings.observer,"DKF"))
                DKF.A = SYS.dsys_obs.A;
                DKF.C = SYS.dsys_obs.C;
                DKF.Bw = SYS.Bw;
                            
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

                % OverWrite derivativeOrder!
                if simulationSettings.batch == true
                    DKF.QuTune = LcurveQus(i);              % Change QuTune
                    DKF.nd = k-1; % Change every loop!
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
                    Uss = c2d(Uss,simulationSettings.dt,mSettings.c2dMethod);
                
                if modelSettings.posMeasurement == true
                    DKF.R = SYS.R*DKF.RTune;                             
                else
                    DKF.R = SYS.R(2:end,2:end)*DKF.RTune;                             
                end

                DKF.uA = Uss.A;
                DKF.uB = Uss.B; 
                DKF.uC = [eye(SYS.nu),zeros(SYS.nu,(DKF.nd)*SYS.nu)];
                
                DKF.D = SYS.dsys_obs.D;
                DKF.B = SYS.dsys_obs.B;
                DKF.Q = eye(size(SYS.Q,1))*DKF.QTune;
                DKF.Qu = eye(SYS.nu)*DKF.QuTune;
        
                if l == 1
                    fprintf(['    Dual Kalman Filter-> \n'...
                     '        QuTune: %1.1e \n'],DKF.QuTune)
                end
            end
            SYS.DKF = DKF;
        
            % GDF setup------------------------------------------------------------
            if any(ismember(simulationSettings.observer,"GDF"))
                GDF.A = SYS.dsys_obs.A;
                GDF.B = SYS.dsys_obs.B;
                GDF.C = SYS.dsys_obs.C;
                GDF.D = SYS.dsys_obs.D;
        
                GDF.Q = eye(size(SYS.Q,1))*GDF.QTune;
                if modelSettings.posMeasurement == true
                    GDF.R = SYS.R*GDF.RTune;
                else
                    GDF.R = SYS.R(2:end,2:end)*GDF.RTune;                             % Snipp off laser measurement
                end

                if l == 1
                    fprintf('    Giljins de Moor filter-> \n')
                end
            end
            SYS.GDF = GDF;
         
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
    %for i = 1:a
    parfor i = 1:a                              
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