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
    
        % MF settings
        MF = [];

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
