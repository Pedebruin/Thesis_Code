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
    modelSettings.patches = [0.11];                             % location of start of patches (y/L)
    modelSettings.nsElementsP = 3;                          % Number of smart elements per patch
    modelSettings.LbElements = 0.1;                         % Preferred length of beam elements (y/L) (will change slightly)
    modelSettings.patchCov = 1e-8;                          % True covariance of patch measurement
    modelSettings.strainRate = false;                        % Measure strain rate instead of strain. 

    %{
        patchCov    low noise: 1e-8
                    high noise: 1e-6
    %}

% Accelerometers
    modelSettings.Acc = [0.95,0.1];                        % Location of accelerometers
    modelSettings.mAcc = 0.005;                              % Mass of accelerometers
    modelSettings.accCov = 0.1;                             % True covariance of accelerometer measurement

    %{
        accCov      low noise: 0.1
                    high noise: 10
    %}
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
    modelSettings.AccError = 1e-3;                      % accelerometer position error

% simulationSettings (Settings for the simulation)%%%%%%%%%%%%%%%%%%%%%%%%%
simulationSettings.simulate = true;                     % Simulate at all?
    simulationSettings.waitBar = false;                      % Give waitbar and cancel option (VERY SLOW)
    simulationSettings.noise = true;                        % Turn on or off all noise!
    simulationSettings.modelError = true;               % Turn on or off al model errors!
    
    % Time settings (Settings for the time and stepsize of the simulation)
    simulationSettings.dt = 1e-3;                       % Sampling time
    simulationSettings.T = 2;                           % Total simulation time
    
    % Input settings (Settings for the input that is used)
    simulationSettings.distInput = 1;                   % Which input is the disturbance?
        simulationSettings.stepTime = [0.5,1];               % Location of input step ([start time], [endtime], [] )
            simulationSettings.stepAmp = 5;             % Step amplitude
        simulationSettings.impulseTime = [];            % Location of input impulse ([time], [])
            simulationSettings.impulseAmp = 10;          % Inpulse amplitude
        simulationSettings.harmonicTime = [];            % Harmonic input start time ([time], [])
            simulationSettings.harmonicFreq = 1;        % Frequency of sinusoidal input ([freq], [])
            simulationSettings.harmonicAmp = 10;         % Frequency input amplitude [Hz]
        simulationSettings.randTime = [];               % random input start time ([time], [])
            simulationSettings.randInt = [-1,1];        % random input interval (uniformly distributed)
    
    % Observer settings (Settings for the observers) 
    simulationSettings.observer = ["AKF"];    % ["MF" "LO" "KF" "AKF" "DKF" "GDF"] Does need to be in order
        simulationSettings.obsOffset = 1e-5;                % Initial state offset (we know its undeformed at the beginning, so probably 0); 
        
        % MF settings
            % The modal filter does not require any settings!

        % LO settings
        LO.poleMovement = 5e-3;                         % 0.5% faster poles

        % KF settings
        KF.stationary = false;	                        % Use stationary KF? 
        KF.QTune = 1e6;
        KF.RTune = 1;

        % AKF settings
        AKF.stationary = false;                         % Use stationary AKF? (DOES NOT WORK YET)
        AKF.derivativeOrder = 0;                        % Higher order derivative? (0:CP, 1:CV, 2:CA)
        AKF.QTune = 1;                                  % Process noise covariance tuning
        AKF.QuTune = 5e3;                               % Input sequence covariance tuning parameter
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
nPatches = length(modelSettings.patches);           % Number of patches
nsElementsP = modelSettings.nsElementsP;            % Number of elements per patch
nAcc = length(modelSettings.Acc);                   % Number of accelerometers

% True model for simulation
SYS_sim = model('True beam model'); % create model object
[SYS_sim, dSYS_sim, modelSettings] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS_sim); % Fill model object 


% Give summary of model
numEl = modelSettings.numEl;
Nmodes = modelSettings.Nmodes;
fprintf(['    # patches: %d \n'...
         '    # accelerometers: %d \n'...
         '    # elements: %d \n'...
         '    # modes: %d \n'],nPatches,nAcc,numEl,Nmodes);

% Observer model with possible modelling error
if simulationSettings.modelError == true
    Beam.rho = Beam.rho*modelSettings.rhoError;                         % Beam mass error
    modelSettings.mAcc = modelSettings.mAcc*modelSettings.mAccError;    % Accelerometer mass error
    modelSettings.Acc = modelSettings.Acc - modelSettings.AccError;     % Acceleromter position error
end

SYS_obs = model('Erroneous beam model');
[SYS_obs, dSYS_obs, modelSettings] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS_obs); % Fill model object

%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulates the system and makes nice plots and stuff. It is not
% really necessary and not the point of this script. But worked wonders when
% debugging the code. 

if simulationSettings.simulate ==  true
    fprintf(['\n Setting up simulation...\n'...
            '    Noise: %d\n'...
            '    Offset: %1.1e\n'...
            '    Model Error: %d\n'],simulationSettings.noise,simulationSettings.obsOffset,simulationSettings.modelError);

    % Snip off laser measurement for observers (As they can't use this)
    dsys_sim = dSYS_sim.sys;                        % Simulated model!
    dsys_obs = dSYS_obs.sys(2:end,1:end);           % Model used for observers!

    % Handy defenitions & Unpacking
    T = simulationSettings.T;
    dt = simulationSettings.dt;
    t = 0:dt:T;

    ny = size(dsys_sim,1);                  % Number of outputs (Total)
    nu = size(dsys_sim,2);                  % Number of inputs
    nq = size(dsys_sim.A,1);                % Number of states

    % Distubrance input generation (function at the bottom)
    Udist = generateInput(simulationSettings);  % Disturbance input

    % Set up waitbar if necessary
    if simulationSettings.waitBar == true
        f = waitbar(0,'1','Name','Running simulation loop...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(f,'canceling',0);
    end

    % MF setup-------------------------------------------------------------
    if any(ismember(simulationSettings.observer,"MF"))
        Psi = pinv(dsys_obs.C);
    end

    % LO setup-------------------------------------------------------------
    if any(ismember(simulationSettings.observer,"LO"))
        poles = eig(dsys_obs.A)*(1-LO.poleMovement);                   % Place poles in Discrete time (TODO)
        LO.L = place(dsys_obs.A',dsys_obs.C',poles)'; 

        fprintf(['    Luenberger Observer-> \n'...
                 '        Pole speed: %3.1f%% increase \n'],LO.poleMovement*100)
    end

    % KF setup-------------------------------------------------------------
    if any(ismember(simulationSettings.observer,"KF"))
        KF.Q = dSYS_obs.Q*KF.QTune;
        KF.R = dSYS_obs.R(2:end,2:end)*KF.RTune;                          % Snip off laser measurement
        KF.S = dSYS_obs.S(:,2:end);
        KF.Bw = dSYS_obs.Bw;
        KF.Dv = dSYS_obs.Dv(2:end,2:end); 
          
        if KF.stationary == true
            % Find stationary kalman gain by solving discrete algebraic ricatti equation
            % (p.162 M.Verhaegen
            [~,Kt,KF.poles] = idare(dsys_obs.A, dsys_obs.C', KF.Q, KF.R, KF.S, []);
            KF.K = Kt';
        else
            % No additional setup for non-stationary kalman filter.
        end

        fprintf(['    Kalman Filter-> \n'...
                 '        Stationary: %d \n'],KF.stationary)
    end

    % AKF setup------------------------------------------------------------
    nd_AKF = AKF.derivativeOrder;
    if any(ismember(simulationSettings.observer,"AKF"))
        % q -> [q;qd;u] ->nq+nu*nd

        AKF.S = [dSYS_obs.S(:,2:end); zeros(nu*(nd_AKF+1),ny-1)];

        % Generate temporary ss model for input dynamics (higher order
        % derivatives can then be modelled!
            uA = [zeros(nu*nd_AKF,nu),eye(nu*nd_AKF);
                    zeros(nu,nu*(nd_AKF+1))];
            uB = [zeros(nd_AKF*nu,nu);
                    eye(nu)];
            uC = zeros(1,nu*nd_AKF);
            uD = [];
            Uss = ss(uA,uB,[],[]);
            Uss = c2d(Uss,dt,modelSettings.c2dMethod);

        AKF.Dv = dSYS_obs.Dv(2:end,2:end);
        AKF.Bw = [dSYS_obs.Bw,zeros(nq,nu);
            zeros(nu*(nd_AKF+1),nq),Uss.B];
        AKF.A = [dsys_obs.A, dsys_obs.B,zeros(nq,nu*nd_AKF);
                zeros(nu*(nd_AKF+1),nq),Uss.A];
        AKF.C = [dsys_obs.C, dsys_obs.D, zeros(ny-1,nu*nd_AKF)];
                                  
        AKF.R = dSYS_obs.R(2:end,2:end)*AKF.RTune;                                       % Do you trust your measurements?
        AKF.Q = dSYS_obs.Q*AKF.QTune;                                           % Do you trust your model?
        AKF.Qu = [1,zeros(1,nPatches);
                  zeros(nPatches,1),zeros(nPatches,nPatches)]*AKF.QuTune;                  % Input covariance!!

        AKF.Q = [dSYS_obs.Q,zeros(nq,nu); % Assemble augmented Q matrix
            zeros(nu,nq),AKF.Qu];
    
        if AKF.stationary == true
            % Find stationary kalman gain by solving discrete algebraic ricatti equation
            % (p.162 M.Verhaegen)
            if nd_AKF > 0 || nu > 1
                error('Ricatti does not work for this augmented system.. (Already in backlog, working on it)')
            end
            [~,Kt,AKF.poles,INFO] = idare(AKF.A, AKF.C', AKF.Q, AKF.R, AKF.S, []);
            AKF.K = Kt';
        else
            % No additional setup for non-stationary augmented kalman filter.
        end

        fprintf(['    Augmented Kalman Filter-> \n'...
                 '        Stationary: %d \n'...
                 '        QuTune: %1.1e \n'],AKF.stationary,AKF.QuTune)
    end

    % DKF setup------------------------------------------------------------
    if any(ismember(simulationSettings.observer,"DKF"))
        nd_DKF = DKF.derivativeOrder;
        DKF.A = dsys_obs.A;
        DKF.B = dsys_obs.B;
        DKF.C = dsys_obs.C;
        DKF.D = dsys_obs.D;

            % Generate temporary ss model for input dynamics (higher order
            % derivatives can then be modelled!
            uA = [zeros(nu*nd_DKF,nu),eye(nu*nd_DKF);
                    zeros(nu,nu*(nd_DKF+1))];
            uB = [zeros(nd_DKF*nu,nu);
                    eye(nu)];
            uC = zeros(1,nu*(nd_DKF+1));
            uC((0:nu-1)*(nd_DKF+1)+1) = ones(1,nu);

            uD = [];
            Uss = ss(uA,uB,uC,[]);
            Uss = c2d(Uss,dt,modelSettings.c2dMethod);
        
        DKF.uA = Uss.A;
        DKF.uB = Uss.B; 
        DKF.uC = Uss.C; % IMPLEMENT THIS IN SIMULATION LOOP! (So does not work yet)

        DKF.Q = dSYS_obs.Q*DKF.QTune;
        DKF.R = dSYS_obs.R(2:end,2:end)*DKF.RTune;                             % Snipp off laser measurement
        DKF.Qu = eye(nu)*DKF.QuTune;

        fprintf(['    Dual Kalman Filter-> \n'...
         '        QuTune: %1.1e \n'],DKF.QuTune)
    end

    % GDF setup------------------------------------------------------------
    if any(ismember(simulationSettings.observer,"GDF"))
        GDF.A = dsys_obs.A;
        GDF.B = dsys_obs.B;
        GDF.C = dsys_obs.C;
        GDF.D = dsys_obs.D;

        GDF.Q = dSYS_obs.Q*GDF.QTune;
        GDF.R = dSYS_obs.R(2:end,2:end)*GDF.RTune;                             % Snipp off laser measurement

        fprintf(['    Giljins de Moor filter-> \n'])
    end
    % Some initialisations!!!
    yfull = zeros(ny,length(t));            % System output
    yfull_MF = zeros(ny,length(t));
    yfull_LO = zeros(ny,length(t)); 
    yfull_KF = zeros(ny,length(t));
    yfull_AKF = zeros(ny,length(t));
    yfull_DKF = zeros(ny,length(t));
    yfull_GDF = zeros(ny,length(t));

    U = zeros(nu,1);                                        % Input vector

    qfull = zeros(nq,length(t));                            % Full state vector over time
    qfull_MF = zeros(nq,length(t));
    qfull_LO = zeros(nq,length(t));                         % Also for the observers
    qfull_KF = zeros(nq,length(t));
    qfull_AKF = zeros(nq+nu*(nd_AKF+1),length(t));          % Augmented state for AKF
    qfull_DKF = zeros(nq,length(t));
    qfull_GDF = zeros(nq,length(t));

    % ufull_AKF can be found in the last two states of qfull_AKF (due tobthe augmented nature)
    ufull_DKF = zeros(nu,length(t));
    ufull_GDF = zeros(nu,length(t));

    q1 = zeros(nq,1);                                       % Normal time simulation    
    q1_LO = ones(nq,1)*simulationSettings.obsOffset;        % Initial state estimate for LO
    q1_KF = ones(nq,1)*simulationSettings.obsOffset;        % Initial state estimate for KF
    q1_AKF = zeros(nq+nu*(nd_AKF+1),1);                         % Initial sate estimate for AKF 
        q1_AKF(1:nq) = ones(nq,1)*simulationSettings.obsOffset; % Only  beam states offset error
    q1_DKF = ones(nq,1)*simulationSettings.obsOffset;       % Initial state estimate for DKF
        u1_DKF = zeros(nu,1);
    q1_GDF = ones(nq,1)*simulationSettings.obsOffset;       % Initial state estimate for GDF
        u1_GDF = zeros(nu,1);

    P1_KF = zeros(nq);                                        % Initial P matrix kalman filter
    P1_AKF = zeros(nq+nu*(nd_AKF+1));
    P1_DKF = zeros(nq);
        Pu1_DKF = zeros(nu);
    Pq1_GDF  = zeros(nq);
        Pu1_GDF = zeros(nu);

    W = mvnrnd(zeros(nq,1),dSYS_obs.Q,T/dt+1)';
    V = mvnrnd(zeros(ny,1),dSYS_obs.R,T/dt+1)';

    fprintf(['    T: %.2f s \n'...
            '    dt: %.4f s\n'...
            '    Time steps: %d \n'],simulationSettings.T, simulationSettings.dt, length(t))
   
    % Simulation loop!!!! Here, the system is propagated.%%%%%%%%%%%%%%%%%%
    fprintf('\n Simulating...\n')
    startTime = tic;
    for i = 1:T/dt+1
        % Shift one time step
        q = q1;             % Previous next state is the current state. (Yes, very deep indeed)
        q_LO = q1_LO;       % Also for the observers
        q_KF = q1_KF;
        q_AKF = q1_AKF;
        q_DKF = q1_DKF;
            u_DKF = u1_DKF;
        q_GDF = q1_GDF;
            u_GDF = u1_GDF;

        P_KF = P1_KF;
        P_AKF = P1_AKF;
        P_DKF = P1_DKF;
            Pu_DKF = Pu1_DKF;
        Pq_GDF = Pq1_GDF;
            Pu_GDF = Pu1_GDF;

        % Simulate system
            % Pick input
            U(simulationSettings.distInput) = Udist(i);

            % Pick noise
            if simulationSettings.noise == true
                w = W(:,i);                         % mvnrnd to allow for different covariances of different sensors
                v = V(:,i);
            else
                w = zeros(nq,1);
                v = zeros(ny,1);
            end
    
            % Propagate discrete time dynamic system
            q1 = dsys_sim.A*q + dsys_sim.B*U + dSYS_obs.Bw*w;           % Propagate dynamics
            y = dsys_sim.C*q + dsys_sim.D*U + dSYS_obs.Dv*v;       % Measurement equation

            % Also save full states for plotting (in q space, so modal)
            qfull(:,i) = q;
            yfull(:,i) = y;

        % Run state estimators
            % MF ----------------------------------------------------------
            if any(ismember(simulationSettings.observer,'MF'))
                % Modal filter
                q_MF = Psi*y(2:end);

                % Save state and output
                qfull_MF(:,i) = q_MF;
                yfull_MF(:,i) = dsys_sim.C*q_MF;
            end

            % LO ----------------------------------------------------------
            if any(ismember(simulationSettings.observer,'LO'))
                q1_LO = dsys_obs.A*q_LO + dsys_obs.B*U + LO.L*(y(2:end)-dsys_obs.C*q_LO-dsys_obs.D*U);
                yfull_LO(:,i) = dsys_sim.C*q_LO + dsys_sim.D*U; % Estimated output
                qfull_LO(:,i) = q_LO;       % Save estimated state
            end 

            % KF ----------------------------------------------------------           
            if any(ismember(simulationSettings.observer,'KF'))
                if KF.stationary == true
                    % Steady state kalman filter (Same as LO, but with kalman gain)
                    q1_KF = dsys_obs.A*q_KF + dsys_obs.B*U + KF.K*(y(2:end)-dsys_obs.C*q_KF-dsys_obs.D*U);

                    % Save state and output
                    qfull_KF(:,i) = q_KF;                    
                    yfull_KF(:,i) = dsys.C*q_KF + dsys.D*U;

                else
                    % Measurement update
                    q_KF = q_KF + P_KF*dsys_obs.C'/(dsys_obs.C*P_KF*dsys_obs.C'+KF.R)*(y(2:end)-dsys_obs.C*q_KF-dsys_obs.D*U);
                    P_KF = P_KF - P_KF*dsys_obs.C'/(dsys_obs.C*P_KF*dsys_obs.C'+KF.R)*dsys_obs.C*P_KF;

                    % Time update
                    q1_KF = dsys_obs.A*q_KF + dsys_obs.B*U;
                    P1_KF = dsys_obs.A*P_KF*dsys_obs.A'+KF.Bw*KF.Q*KF.Bw';

                    % Save state and output
                    qfull_KF(:,i) = q_KF;                    
                    yfull_KF(:,i) = dsys_sim.C*q_KF + dsys_sim.D*U;

                end
            end

            % AKF----------------------------------------------------------
            if any(ismember(simulationSettings.observer,'AKF'))
                if AKF.stationary == true
                    % Steady state kalman filter (Same as LO, but with kalman gain)
                    q1_AKF = AKF.A*q_AKF + AKF.K*(y(2:end,i)-AKF.C*q_AKF);

                    % Save state and output
                    qfull_AKF(:,i) = q_AKF;                    
                    yfull_AKF(:,i) = [dsys.C(1,:),zeros(1,nu);
                                    AKF.C]*q_AKF;                   
                else
                    % Measurement update
                    q_AKF = q_AKF + P_AKF*AKF.C'/(AKF.C*P_AKF*AKF.C'+AKF.R)*(y(2:end)-AKF.C*q_AKF);
                    P_AKF = P_AKF - P_AKF*AKF.C'/(AKF.C*P_AKF*AKF.C'+AKF.R)*AKF.C*P_AKF;

                    % Time update
                    q1_AKF = AKF.A*q_AKF;
                    P1_AKF = AKF.A*P_AKF*AKF.A'+AKF.Bw*AKF.Q*AKF.Bw';

                    % Save state and output
                    qfull_AKF(:,i) = q_AKF;                    
                    yfull_AKF(:,i) = [dsys_sim.C(1,:),zeros(1,nu*(nd_AKF+1));
                                    AKF.C]*q_AKF;
                end
            end

            % DKF----------------------------------------------------------
            if any(ismember(simulationSettings.observer,'DKF'))
                % Measurement update of INPUT estimate
                Ku_DKF = Pu_DKF*DKF.D'/(DKF.D*Pu_DKF*DKF.D' + DKF.R);

                u_DKF = u_DKF + Ku_DKF*(y(2:end)-DKF.C*q_DKF-DKF.D*u_DKF);
                Pu_DKF = Pu_DKF - Ku_DKF*DKF.D*Pu_DKF;

                % Measurement update STATE estimate
                K_DKF = P_DKF*DKF.C'/(DKF.C*P_DKF*DKF.C' + DKF.R);

                q_DKF = q_DKF + K_DKF*(y(2:end)-DKF.C*q_DKF-DKF.D*u_DKF);
                P_DKF = P_DKF - K_DKF*DKF.C*P_DKF;

                % Time update of INPUT estimate
                u1_DKF = u_DKF;                 % Random walk!
                Pu1_DKF = Pu_DKF + DKF.Qu;

                % Time update STATE estimate
                q1_DKF = DKF.A*q_DKF + DKF.B*u_DKF;  % Dynamics propagation using estimated input
                P1_DKF = DKF.A'*P_DKF*DKF.A + DKF.Q;

                % Save state,estimated input and output 
                qfull_DKF(:,i) = q_DKF;    
                ufull_DKF(:,i) = u_DKF;
                yfull_DKF(:,i) = [dsys_sim.C(1,:);
                                DKF.C]*q_DKF; 
               
            end
            
            % GDF----------------------------------------------------------
            if any(ismember(simulationSettings.observer,'GDF'))
                % Input estimation
                Rt_GDF = GDF.C*Pq_GDF*GDF.C' + GDF.R;
                M_GDF = (GDF.D'/Rt_GDF*GDF.D)\GDF.D'/Rt_GDF;
                u_GDF = M_GDF*(y(2:end)-GDF.C*q_GDF);
                Pu_GDF = eye(nu)/(GDF.D'/Rt_GDF*GDF.D);

                % Measurement update
                K_GDF = Pq_GDF*GDF.C'/GDF.R;
                q_GDF = q_GDF + K_GDF*(y(2:end)-GDF.C*q_GDF-GDF.D*u_GDF);
                Pq_GDF = Pq_GDF - K_GDF*(Rt_GDF - GDF.D*Pu_GDF*GDF.D')*K_GDF';
                Pqu_GDF = -K_GDF*GDF.D*Pu_GDF;

                % Time update
                q1_GDF = GDF.A*q_GDF + GDF.B*u_GDF;
                Pq1_GDF = [GDF.A GDF.B]*[Pq_GDF ,Pqu_GDF; Pqu_GDF' Pu_GDF]*[GDF.A'; GDF.B'] + GDF.Q;

                % Save state, estimated input and output
                qfull_GDF(:,i) = q_GDF;    
                ufull_GDF(:,i) = u_GDF;
                yfull_GDF(:,i) = [dsys_sim.C(1,:);
                                GDF.C]*q_GDF; 
            end

        % Update waitbar and check cancel button
        if simulationSettings.waitBar == true
            waitbar(i/(T/dt+1),f,sprintf('Step %d/%d',i,T/dt+1))
            if getappdata(f,'canceling')
                break
            end
        end
    end
    elapsed = toc(startTime);
    fprintf('    Simulation time: %.2f s \n',elapsed)

    if simulationSettings.waitBar == true
        delete(f);
    end

    % Save simulation data in the larger model struct. 
    dSYS_sim.simulationData.t = t;
    dSYS_sim.simulationData.Udist = Udist;
    
    dSYS_sim.simulationData.qfull = qfull;
    dSYS_sim.simulationData.qfull_MF = qfull_MF;
    dSYS_sim.simulationData.qfull_LO = qfull_LO;
    dSYS_sim.simulationData.qfull_KF = qfull_KF;
    dSYS_sim.simulationData.qfull_AKF = qfull_AKF;
    dSYS_sim.simulationData.qfull_DKF = qfull_DKF;
    dSYS_sim.simulationData.qfull_GDF = qfull_GDF;

    dSYS_sim.simulationData.yfull = yfull;
    dSYS_sim.simulationData.yfull_MF = yfull_MF;
    dSYS_sim.simulationData.yfull_LO = yfull_LO;
    dSYS_sim.simulationData.yfull_KF = yfull_KF;
    dSYS_sim.simulationData.yfull_AKF = yfull_AKF;
    dSYS_sim.simulationData.yfull_DKF = yfull_DKF;
    dSYS_sim.simulationData.yfull_GDF = yfull_GDF;

    dSYS_sim.simulationData.ufull_DKF = ufull_DKF;
    dSYS_sim.simulationData.ufull_GDF = ufull_GDF;

    % Make a nice plot of the simulation!
    if plotSettings.plot == true
        fprintf('\n Plotting simulation... \n')
        simPlots = plotter(dSYS_sim);
    end
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