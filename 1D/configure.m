function [Beam, sBeam, modelSettings,...
        simulationSettings, plotSettings, simulinkSettings,...
        LO,KF,AKF,DKF,GDF,...
        Cd, configuredExternally] = configure(External,varargin)
    %{
        This file serves as the configuration for both the matlab and simulink
        files. All settings are placed here!
    %}
    
    %% Configure externally inputted variables
    
    %% Model Parameters & Settings
    % To evaluate different beam properties

    % General parameters that can be optimised!
    L = 99e-3;
    b = 12.7e-3;
    h = 0.4e-3;
    E = 200E9;
    rho = 7850;

    stageMass = 100e-3;
    patchL = 22.5e-3;
    patches = [2e-3/L];
    simulate = true;

    QTuneAKF= 1.86e-9;
    RTuneAKF= 2.26e4;
    QuTuned0AKF= 1e10;

    QTuneDKF = 2.2723e-13;
    RTuneDKF = 1.4073e4;
    QuTuned0DKF = 1e1; 

    QTuneGDF = 2.2375e-5;
    RTuneGDF = 2.2885e14; 
    
    % Overwrite when configured externally
    if External == true && length(varargin) >= 1 % Put default values in
        parameters = varargin{1};
        values = varargin{2};
        baseMaterial = varargin{3};
        switch baseMaterial
            case "Ssteel"
%                 (Commented out to make sure no values are
%                 overwritten incorrectly. In this case, the parameters at
%                 the top of the file are used for SS. 
%                 L = 99e-3; 
%                 b = 12.7e-3;
%                 h = 0.4e-3;
%                 E = 200E9;
%                 rho = 7850;
            case "Aluminium"
                L = 98.3e-3;
                b = 15e-3;
                h = 1e-3;
                E = 68.5e9;
                rho = 2700;
            otherwise
                error('No base material selected')
        end
        
        for i = 1:length(parameters)            % Change values as required
            parameter = parameters(i);
            switch parameter
                case "L" % Length
                    L = values(i);
                case "b" % Width
                    b = values(i);
                case "h" % Thickness
                    h = values(i);
                case "E" % E modulus
                    E = values(i);
                case "rho" % Density
                    rho = values(i);
                case "stageMass" % stageMass
                    stageMass = values(i);
                case "patchL" % Patch length
                    patchL = values(i);
                case "patches" % Patch location
                    patches = values(i);
                case "simulate" % Simulate or not?
                    simulate = values(i); 
                case "QTuneAKF" % Covariance tuning parameter AKF
                    QTuneAKF = values(i);
                case "RTuneAKF" % 
                    RTuneAKF = values(i);
                case "QuTuned0AKF"
                    QuTuned0AKF = values(i);
                case "QTuneDKF"
                    QTuneDKF = values(i);
                case "RTuneDKF"
                    RTuneDKF = values(i);
                case "QuTuned0DKF"
                    QuTuned0DKF = values(i);
                case "QTuneGDF"
                    QTuneGDF = values(i);
                case "RTuneGDF"
                    RTuneGDF = values(i);
            end  
        end
    else

    end
    
    % Model settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % beam element parameters (This is a beam element)
        Beam = element('Beam');
        Beam.h = h;                         % m thickness of beam
        Beam.b = b;                         % m width of beam
        Beam.E = E;                      % E modulus of beam
        Beam.mu = 0.334;                    % Poisson's ratio of beam
        Beam.rho = rho;                    % Mass density of beam Kg/M^3
        Beam.zeta = 0.002;                   % Modal damping coefficient of beam 
    
    % Stage mass (For 1 flexure, Only for testing purposes)
        modelSettings.stageMass = stageMass;
    
    % Smart patches (Piezo or strain)
        modelSettings.patches = patches;                                 % location of start of patches (y/L)
        modelSettings.gauges = "resistive";                         % "resistive" or "piezo" gauges? 
        modelSettings.nsElementsP = 3;                              % Number of smart elements per patch
        modelSettings.strainRate = false;                            % Measure strain rate instead of strain. (For current amplifiers)
        modelSettings.posMeasurement = false;                        % Use position measurement for the observers?
    
        sBeam = copy(Beam); sBeam.name = 'sBeam';   
        switch modelSettings.gauges
            case {'piezo'}                                          % When using piezo patches
                sBeam.ph = 0.5e-3;                                  % m 
                sBeam.pb = b;                                       % m
                sBeam.pE = 5.2e10;                      
                sBeam.pmu = 0.334;
                sBeam.prho = 7800;
                modelSettings.patchCov = 1e-1;                    	% True covariance of patch measurement (log10)
            case {'resistive'}                                      % When using resistive strain gauges
                sBeam.ph = 0.025e-3;                                % m 
                sBeam.pb = b;                                       % m
                sBeam.pE = 10;                      
                sBeam.pmu = 0.334;
                sBeam.prho = 1000;
                modelSettings.patchCov = 1.095e-16;                       % True covariance of patch measurement (log10)
                % Value gets overwritten with actual covariance in MAIN
                % line 154 when using real measurements. 
        end
        sBeam.L = patchL/modelSettings.nsElementsP;                 % Length of smart beam element
    
    % Accelerometers
        modelSettings.Acc = [1];                            	    % Location of accelerometers
        modelSettings.mAcc = 0;%0.25e-3;                               % Mass of accelerometers
        modelSettings.accCov = 0.0587;                                   % True covariance of accelerometer measurement
        % Value gets overwritten with actual covariance in MAIN
        % line 160 when using true real measurements.

    % General Modelling 
        modelSettings.LbElements = 0.1;                             % Preferred length of beam elements (y/L) (will change slightly)
        modelSettings.Nmodes = 2;                                   % Number of modes to be modelled (Can't be larger then the amount of nodes)
        modelSettings.measurementHeight = 1;                        % Height of the measurement
        modelSettings.forceHeight = 1;                              % Height of the input force.
        
        modelSettings.wcov = 0;                              % Process noise covariance
        modelSettings.laserCov = 0;%1.4210e-12;                          % Laser covariance                                         
        % Value gets overwritten with actual covariance in MAIN
        % line 162 when using real measurements. 
        
        modelSettings.L = L;                                    
        modelSettings.b = b;
        
        modelSettings.rhoError = 1.05;                      % []% density error
        modelSettings.mAccError = 1.05;                     % []% accelerometer mass error
        modelSettings.AccError = 1e-3;                      % accelerometer position error
        
        modelSettings.c2dMethod = 'ZOH';
        
        %% Observer settings
        simulationSettings.observer = ["MF" "LO" "KF" "AKF" "DKF" "GDF"];    % ["MF" "LO" "KF" "AKF" "DKF" "GDF"] Does need to be in order
            simulationSettings.obsOffset = 1e-10;               % Initial state offset (we know its undeformed at the beginning, so probably 0); 
        
            % MF settings
            MF = [];
    
            % LO settings
            LO.poleMovement = 5e-3;                         % 0.5% faster poles
    
            % KF settings
            KF.stationary = true;	                        % Use stationary KF? 
            KF.QTune = 1.86e-9;
            KF.RTune = 2.26e4;
    
            % AKF settings
            AKF.stationary = false;                         % Use stationary AKF? (DOES NOT WORK YET)
            AKF.derivativeOrder = 0;                        % % Derivative order when not in batch mode (0:CP, 1:CV, 2:CA)
            AKF.QTune = QTuneAKF; % eye*QTune                      % Process noise covariance tuning
            AKF.RTune = RTuneAKF; % R*RTune                        % Measurement noise covariance tuning
            AKF.P0 = 1e3;
            AKF.Pu0 = 1e3;
    
            AKF.QuTuned0 = QuTuned0AKF;                             % For the zeroth derivative 
            AKF.QuTuned1 = [];                             % For the first derivative
            AKF.QuTuned2 = [];                             % For the second derivative
    
            % DKF settings
            DKF.derivativeOrder = 0;                                                                            % Derivative order when not in batch mode
            DKF.QTune = QTuneDKF; % eye*QTune
            DKF.RTune = RTuneDKF; % R*RTune
            DKF.P0 = 1e3;
            DKF.Pu0 = 1e3;
    
            DKF.QuTuned0 = QuTuned0DKF;                                                                         
            DKF.QuTuned1 = [];
            DKF.QuTuned2 = [];
            
            % GDF settings
            GDF.QTune = QTuneGDF; % eye*QTune
            GDF.RTune = RTuneGDF;  
            GDF.Pq0 = 1e3;
            GDF.Pu0 = 1e3;
            GDF.Pqu0 = 1e3;
            
    %% Matlab Simulation settings 
    simulationSettings.simulate = simulate;                     % Simulate at all?
              
        simulationSettings.waitBar = false;                      % Give waitbar and cancel option (VERY SLOW)
                
        simulationSettings.data = "Real";                       % Use Real or Simulated data (Real/Simulated)
            simulationSettings.dataset = 'FINALOL2';   % Long version
            
            % Use the true or faked accelerometer measurements?
            simulationSettings.trueAccMeasurement = true;      % Use true or faked accelerometer settings (only for real dataset)
            
            % Turn on and off all errors!
            simulationSettings.noise = true;                   % Turn on or off all noise!
            simulationSettings.modelError = true;              % Turn on or off al model errors!
            simulationSettings.initialError = true;            % Turn on or off all initial state errors!
           
            % Time settings (Settings for the time and stepsize of the simulation)
            simulationSettings.dt = 1e-4;                       % Sampling time
            simulationSettings.T = 10;                           % Total simulation time
            
            % Input settings (Settings for the input that is used)
            simulationSettings.distInput = true;                   % Which input is the disturbance?
                simulationSettings.stepTime = [0.5,1.5];               % Location of input step ([start time], [endtime], [] )
                    simulationSettings.stepAmp = 0.05;             % Step amplitude
                simulationSettings.impulseTime = [];            % Location of input impulse ([time], [])
                    simulationSettings.impulseAmp = 1;          % Inpulse amplitude
                simulationSettings.harmonicTime = [];            % Harmonic input start time ([time], [])
                    simulationSettings.harmonicFreq = 538;        % Frequency of sinusoidal input ([freq], [])
                    simulationSettings.harmonicAmp = 1;         % Frequency input amplitude [Hz]
                simulationSettings.randTime = [];               % random input start time ([time], [])
                    simulationSettings.randInt = [-1,1]*1e-2;        % random input interval (uniformly distributed)
            simulationSettings.lowpassInput = false;
                simulationSettings.cutoffFrequency = 150;       % Hz

            % Tune filters on real dataset to finetune parameters (Only
            % works with real dataset, might take a while. If not used,
            % standard parameters are used. 
            simulationSettings.tuneFilters = false;
                simulationSettings.maxFunEvals = 10000;
                simulationSettings.filtersToTune = ["DKF"];

            simulationSettings.saveResult = false;

        % Batch run settings
        %{
            Batch mode runs the script multiple times to analyse meta
            results like the influence of different tuning parameters, or
            noise on performance. The script is currently configured to
            vary QuTune to generate an L-curve plot. 
        %}
        simulationSettings.batch = false;                   % Run simulation in batch?? (A lot of times)
            simulationSettings.nDerivatives = 0;            % Highest derivative order? (Also does all lower derivatives) 
            simulationSettings.QuMin = 1e-15;                 % Loop from
            simulationSettings.QuMax = 1e15;                % Loop to
            simulationSettings.nLcurve = 30;                 % In .. steps
    
        % Run each simulation multiple times for noise realisations?
        simulationSettings.monteCarlo = false;           % Run every simulation multiple times to get a mean and monteCarlo interval
            simulationSettings.iterations = 5;          % Amount of times every simulation is ran. 
            % Also change for to parfor in MAIN line 321!! (Can't be done
            % automatically)

        % Parallel computing settings
        simulationSettings.parallel = false;                % Run simulations in parallel?
            simulationSettings.nWorkers = 4;                % Number of cores to run this on?
    
    %% Matlab plot settings
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
            plotSettings.states = 2;                            % First # states to be plotted
    
        plotSettings.modes = 0;                                 % Plot the mode shapes?? [number of modes]
            plotSettings.modeAmp = 5e-4;                        % Amplification factor for plot. 
    
        plotSettings.alphaMin = 0.3;                            % start of alpha interp

        plotSettings.plotBodePeak = false;                      % Some other plot
    
    %% Simulink parameters & settings
    % General
    simulinkSettings.measurementSampleTime = 1e-4;
    simulinkSettings.PIDSampleTime = 1e-4;
    simulinkSettings.laserResolution = 39.5e-9; %Step size of the laser encoder
    simulinkSettings.maxEncoderPosition = 0.3e5;
    
    % Strain settings
    simulinkSettings.strainGain = -1.3636/(1e6);%-1/(6.0143e5);
    modelSettings.strainCorrectionGain = 0.68;
    
    % Accelerometer settings
    simulinkSettings.accSensitivity = -80e-3/9.813*((2^16)*3.6)/22.95;  %V/(m/s^2)
    simulationSettings.accSensitivity = simulinkSettings.accSensitivity;
    
    % Identification
    simulinkSettings.matchingGain = 1/db2mag(-48.6);
    simulinkSettings.inputGain = 100;
    simulinkSettings.chirpInitFreq = 100;
    simulinkSettings.chirpTargetFreq = 300;
    simulinkSettings.chirpTargetTime = 60;
    
    % Controller settings
    simulinkSettings.controlBandwidth = 1e3;%750;
    
    % Input generation
    simulinkSettings.referenceSampleTime = 1e-4;
    simulinkSettings.referenceAmplitude = 3e-4; %  m Setpoint position
    simulinkSettings.referenceFrequency = 10;%538; % Hz

    simulinkSettings.impulseTime = 10;
    simulinkSettings.impulseSamples = 10;

    % When constant force is given
    simulinkSettings.inputForce = 0.05; %N per flexure
    
    % Get controller from file
    Cd = load('Cd');
    Cd = Cd.Cd;

    % Tell wether the configuration is called externally (outside of MAIN)
    configuredExternally = External;
end