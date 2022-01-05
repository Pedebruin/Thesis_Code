%{
This is a 1D Euler bernoulli simulation with piëzo patches

The beam model is generated based on Euler Bernoulli beam theory. The
piezoelec coupling is implemented as in the paper from 

Aktas, K. G.
Esen, I.
Doi:10.48084/etasr.3949

The author of this code is Pim de Bruin. 

The main model parameters are found at the top of the script. The feedback
and piezo parameters are found around line 260 and 200 respectively. 

The patches can be placed through modelSettings.patches. This is a vector
with the bottom location of each patch in natural coordinates. the length
of this vector determines the amount of patches. The size of the patches is 
determined through patchL. 

The rest of the parameters are relatively self explanatory. If you have any
questions or find an error (very probable) dont hesitate to send me a
message!

XXX Pim

This script requires:
    Symbolic math toolbox 
    Control system toolbox (For feedback command)
    Robust Control toolbox (For Hinf command)
%}

clear
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

disp('Setting up...')
%% Parameters & Settings
% Basic settings
L = 370e-3;     % Total length
b = 40e-3;      % Total width
patchL = 50e-3; % patch length

% Model settings
% Smart patches (Piezo)
modelSettings.patches = [0,0.5];                             % location of start of patches (y/L)
modelSettings.nsElementsP = 3;                          % Number of smart elements per patch
modelSettings.LbElements = 0.1;                         % Preferred length of beam elements (y/L) (will change slightly)

% Accelerometers
modelSettings.Acc = [0.45,0.85];                            % Location of accelerometers
modelSettings.mAcc = 0.01;                              % Mass of accelerometers

% Strain patches (TODO)

% Modelling 
modelSettings.Nmodes = 10;                               % Number of modes to be modelled (Can't be larger then the amount of nodes)
modelSettings.measurementHeight = 1;                 % Height of the measurement
modelSettings.forceHeight = 0.5;                    % Height of the input force. 

modelSettings.L = L;                                    
modelSettings.b = b;

% plotSettings
plotSettings.plot = true;                         % Plot the simulation?
    plotSettings.plotNodes = true;                          % Plot the nodes
        plotSettings.nodeNumbers = false;                   % Plot the node numbers
    plotSettings.elementNumbers = true;                     % Plot the element numbers
    plotSettings.sensor = true;                             % Plot the sensor in beam plot (red line)
    plotSettings.Input = false;                              % Plot given force input in the sensor plot
    plotSettings.accelerometers = true;                     % Plot accelerometers?
        plotSettings.accNumbers = true;            % Plot acceleromter numbers?

% simulationSettings
simulationSettings.simulate = true;                     % Simulate at all?
simulationSettings.waitBar = false;                      % Give waitbar and cancel (takes longer, but good in debugging)
    % Time settings
    simulationSettings.dt = 1e-3;                       % Sampling time
    simulationSettings.T = 2;                           % Total simulation time
    
    % Input settings
    simulationSettings.distInput = 1;                   % Which input is the disturbance?
        simulationSettings.stepTime = [];          % Location of input step ([time], [time endtime], [] )
            simulationSettings.stepAmp = 1;             % Step amplitude
        simulationSettings.impulseTime = [];           % Location of input impulse ([time], [])
            simulationSettings.impulseAmp = 1;          % Inpulse amplitude
        simulationSettings.harmonicTime = [0.1];            % Harmonic input start time ([time], [])
            simulationSettings.harmonicFreq = 1;        % Frequency of sinusoidal input ([freq], [])
            simulationSettings.harmonicAmp = 1;         % Frequency input amplitude [Hz]
        simulationSettings.randTime = [];                % random input start time ([time], [])
            simulationSettings.randInt = [-1,1];        % random input interval (uniformly distributed)


% beam element parameters (for every beam element)
Beam = element('Beam');
Beam.h = 1e-3;                      % m thickness of beam
Beam.b = b;                         % m width of beam
Beam.E = 70E9;                      % E modulus of beam
Beam.mu = 0.334;                    % Poisson's ratio of beam
Beam.rho = 2710;                    % Mass density of beam
Beam.zeta = 0.01;                   % Modal damping coefficient of beam

% Smart beam parameters :PI P-876 DuraAct (for every smart element)
sBeam = copy(Beam); sBeam.name = 'sBeam';       
sBeam.L = patchL/modelSettings.nsElementsP;     % Length of smart beam element
sBeam.ph = 0.5e-3;                              % m 
sBeam.pb = b;                                   % m
sBeam.pE = 5.2e10;                      
sBeam.pmu = 0.334;
sBeam.prho = 7800;

% strain element parameters (TODO)

%% Matrix setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPatches = length(modelSettings.patches);           % Number of patches
nsElementsP = modelSettings.nsElementsP;            % Number of elements per patch
nAcc = length(modelSettings.Acc);                   % Number of accelerometers

% Get element lengths! (Function at the bottom)
[Ls,sElements] = getLengths(modelSettings,sBeam);   % Vector with for every element its length, and which elements are smart
modelSettings.sElements = find(sElements);          % Which elements are smart?

% Find accelerometer locations
accHeights = L*modelSettings.Acc;
nodeHeights = [0,cumsum(Ls)];
accPos = zeros(nAcc,2);         % accelerometer positions [Element, eta]

for i = 1:nAcc
    accHeight = accHeights(i);
    diff = nodeHeights-accHeight;
    upperNode = find(round(diff,4)>=0,1);           % below this node
    lowerNode = upperNode-1;                        % Above this node
    accElement = lowerNode;                         % This element!

    % Natural element coordinate (How far along this element, for integration)
    eta = (accHeight-nodeHeights(lowerNode))/(nodeHeights(upperNode)-nodeHeights(lowerNode));
    accPos(i,1) = accElement;
    accPos(i,2)= eta;
    accPos(i,3)= i;
end
modelSettings.accElements = accPos;                 % Where are the accelerometers placed?

% Assembly
numEl = length(Ls); modelSettings.numEl = numEl;
numNodes = numEl+1; modelSettings.numNodes = numEl+1;
Nmodes = modelSettings.Nmodes;

K = zeros(numNodes*2);                      
M = zeros(numNodes*2);
elements = element.empty(numEl,0);          % Empty element array going to store all elements. 

% Loops over all elements, constructs their elemental K and M matrices and
% stores them in element objects aswell (class defenition in element.m).
% Also constructs the system K and M matrices by putting elemental matrices
% in there. 
for i = 1:numEl
    % Elemental coordinates: qe = [v1,theta1,v2,theta2]';
    n1 = i;     % Starting node
    n2 = i+1;   % Ending node
  
    if sum(ismember(modelSettings.sElements,i)) > 0 % Is smart element?
        El = copy(sBeam);                           % Inherit all smart element parameters
        El.sBeam = true;                            % Tell it its a smart element
    else                                            % Or normal beam element?
        El = copy(Beam);                            % Inherit regular beam element parameters
    end
    
        El.L = Ls(i);   % Give element correct length
        El.ainv();      % Calculate Ainv for interpolation (speed)
        El.number = i;  % Give element correct element number
        
        % Calculate stiffness (for normal beam ph & pb will be zero, so no additional terms). 
        Ib = El.b*El.h^3/12;                        % Beam only
        Ip = El.pb*El.ph^3/12 + El.pb*El.ph*(El.ph+El.h)^2/4; % Patch only
        
        EI = El.E*Ib+2*El.pE*Ip;                     % Equivalent EI and Rho*A (for normal element, extra term is 0)
        rhoA = El.b*(El.rho*El.h+2*El.prho*El.ph);

        % Elemental K and M matrices
        Ke = EI/El.L^3*    [12, 6*El.L, -12, 6*El.L;
                            6*El.L, 4*El.L^2, -6*El.L, 2*El.L^2;
                            -12, -6*El.L, 12, -6*El.L;
                            6*El.L, 2*El.L^2, -6*El.L, 4*El.L^2];

        Me = rhoA*El.L/420* [156, 22*El.L, 54, -13*El.L;
                            22*El.L, 4*El.L^2, 13*El.L, -3*El.L^2;
                            54, 13*El.L, 156, -22*El.L;
                            -13*El.L, -3*El.L^2, -22*El.L, 4*El.L^2];   

    % Add point masses for accelerometers
    if ismember(El.number,accPos(:,1))              % Does this element carry an accelerometer?
        accNumber = find(accPos(:,1) == El.number); % Number of current accelerometer
        eta = accPos(accNumber,2);

        alpha = eta*El.L;                                   % Distance along the element
        for j = 1:length(alpha)                             % Possibly multiple accelerometers per element!
            phi = [1, alpha(j), alpha(j)^2, alpha(j)^3]*El.Ainv;         % Shape functions at location
            Meacc = modelSettings.mAcc*(phi'*phi);              % Mass matrix due to accelerometer
            Me = Me + Meacc;                                    % Add to elemental mass matrix! (might need to check)
        end

        % Also tell the element it has (an) accelerometer(s)
        El.Acc = true;
        El.accEta = eta;
        El.accNumber = accNumber;
    end
    
    % Apply boundary condition to first node
    if n1 == 1                                      % If the node we are looking at is the first one
        % Prescribe 0 to displacement and rotation
        Ke(1:2,:) = zeros(2,4);                   
        Ke(:,1:2) = zeros(4,2);
        Ke(1:2,1:2) = eye(2);                   % Fix both rotation and displacement
        
        % Prescribe 0 to acceleration aswell
        Me(1:2,:) = zeros(2,4);
        Me(:,1:2) = zeros(4,2);
        Me(1:2,1:2) = eye(2);
    end

    % Assembly. Since everything is nicely connected in sequence (1to2to3..) , the
    % elemenental matrices stay nicely together in a diagonal fasion.
    K(n1*2-1:n2*2,n1*2-1:n2*2) = K(n1*2-1:n2*2,n1*2-1:n2*2) + Ke;   % K
    M(n1*2-1:n2*2,n1*2-1:n2*2) = M(n1*2-1:n2*2,n1*2-1:n2*2) + Me;   % M
    
    % Put information in element objects aswell (mainly for plotting, also for structuring)
    x1 = 0;                 % x position first node of element
    x2 = 0;                 % x position second node of element
    y1 = sum(Ls(1:i-1));    % y position first node of element
    y2 = sum(Ls(1:i));      % y position second node of element
    t1 = 0;                 % theta of first node of element
    t2 = 0;                 % theta of second node of element
    initPos = [x1, x2;
                y1, y2;
                t1, t2];
    
    elements(i) = El.assign(i,[n1,n2],initPos,Ke,Me);       % Write everything to the element
end


% Find mode shapes and natural frequencies
if modelSettings.Nmodes <= numEl
        [Phi,omega2] = eig(K,M);
else
    error('Number of elements is smaller then the number of modes requested. Modal decomposition will not work!')
end

Phi = Phi(:,1:Nmodes+2);                % Select only Nmodes+2 modes
omega2 = omega2(1:Nmodes+2,1:Nmodes+2); % Select only Nmodes+2 natural frequencies

%Snip off weird modes 
Phi = Phi(:,3:end);                     % First 2 modes are not correct, so snip them off
omega2 = omega2(3:end,3:end);           % Same for natural frequencies

% Check first natural frequency
f = diag(sqrt(omega2)/(2*pi));
I = (Beam.h^3)/12;
S = Beam.h;         
lambda = 1.87510407;
analyticalf1 = lambda^2/(2*pi*L^2)*sqrt(Beam.E*I/(Beam.rho*S));

delta = 1*L^3/(3*Beam.E*I);
  
if abs(f(1) - analyticalf1) > 0.001 && nPatches == 0 && nAcc == 0 % Check only when there are no patches
    % Is only checked if there are no patches or accelerometers present!
    warning('Eigenfrequencies do not correspond')
end


%% Dynamical matrices setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi = Phi/(Phi'*M*Phi);     % Normalise w.r.t. mass matrix

% Piezo in & outputs
    % Sensor equation according to K. Aktas
    [intS,~] = shapeFunctions();
    H = 1e6;                                        % Arbitrary gain (needs setup validation)
    z = sBeam.h/2+sBeam.ph;                         % Effective height
    d31 = -180e-12;                                 % Piezo coupling d
    s11 = 16.1e-12;
    e31 = d31/s11;                                  % Piezo coupling e
    w = sBeam.pb;                                   % width of beam
    S = H*z*e31*w*intS;                       % Need to check this again (derivation)!
    
    % Actuator equation according to K. Aktas
    [~,intG] = shapeFunctions();
    Ep = sBeam.pE;
    d31 =  -180e-12;
    w = sBeam.pb;
    zbar = (sBeam.ph+sBeam.h)/2;
    %intG = [0,-1,0,1];
    G = Ep*d31*w*zbar*intG';                  % While we're at it, check this derivation aswell!
    
    % Voltage input B matrix (in d coordinates) and output matrix C for the
    % piezo patches!
    nsElements = length(modelSettings.sElements);
    Bg = zeros(numNodes*2,nPatches);
    Cs = zeros(nPatches,numNodes*2);
    for i = 1:nPatches
        for j = 1:nsElementsP
            el = modelSettings.sElements((i-1)*nsElementsP + j);
            n1 = elements(el).neighbours(1);
    
            Bg(n1*2-1:n1*2+2,i) = Bg(n1*2-1:n1*2+2,i) + G;
            Cs(i,n1*2-1:n1*2+2) = Cs(i,n1*2-1:n1*2+2) + S;
        end
    end

% Bmatrices
    % External input B matrix (in d coordinates)
        % Find interpolate nodes
%         height = modelSettings.measurementHeight*L;
%         pos = unique([elements.pos]');
%         diff = pos-height;
%         upperNode = find(round(diff,4)>=0,1);
%         interpNodes = [upperNode-1,upperNode];  
%         interpEl = upperNode-1;
% 
%         % local coordinate alpha
%         alpha = height - pos(interpNodes(1));                       % alpha = eta*obj.L
%         N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
%         
%         Bext(interpNodes(1)*2-1:interpNodes(1)*2+2) = N;
        Bext = zeros(numNodes*2,1);
        fNode = numNodes;                                      % force Node
        Bext(fNode*2-1) = 1;                                % Force at last node in x direction 

% C matrices    
    % Measurement C matrix for laser measurement
        % Find interpolate nodes
        height = modelSettings.measurementHeight*L;
        pos = unique([elements.pos]');
        diff = pos-height;
        upperNode = find(round(diff,4)>=0,1);
        interpNodes = [upperNode-1,upperNode];  
        interpEl = upperNode-1;
        
        % local coordinate alpha
        alpha = height - pos(interpNodes(1));                       % alpha = eta*obj.L
        N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
        Cmeasurement = zeros(1,numNodes*2);                         
        Cmeasurement(interpNodes(1)*2-1:interpNodes(1)*2+2) = N;    % Construct C from this (cubic interpolation)
    
    % Acceleration measurement C matrix
        Cacc = zeros(nAcc,numNodes*2);
        for i = 1:nAcc
            interpEl = accPos(i,1);
            alpha = accPos(i,2)*elements(interpEl).L;
            interpNodes = [interpEl,interpEl+1];
            N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
            Cacc(i,interpNodes(1)*2-1:interpNodes(1)*2+2) = N;
        end
       
% Damping matrix 
    Cmodal = 2*Beam.zeta*sqrt(omega2);

%% Construct full system 
A = [zeros(Nmodes),eye(Nmodes);
    -omega2,-Cmodal];
B = [zeros(Nmodes,nPatches+1);                          
    Phi'*[Bext,Bg]];                             % [Force input, Piezo inputs]
C = [Cmeasurement*Phi, zeros(1,Nmodes);                 % Laser measurement
    Cs*Phi,zeros(nPatches,Nmodes);                      % Piezo outputs
    -Cacc*Phi*omega2, -Cacc*Phi*Cmodal];                % Accelerometer outputs
D = [zeros(nPatches+1,size(B,2));                       % No throughput laser,base and piezo measurements
    Cacc*(Phi*Phi')*[Bext,Bg]];                  % Throughput for acceleration measurements!

sys = ss(A,B,C,D);

% Put the system together with other information in a larger model object. 
SYS = model('Continuous time beam model'); % create model object called SYS!
    SYS.elements = elements;
    SYS.analyticalf1 = analyticalf1;
    SYS.numEl = numEl;
    SYS.numNodes = numNodes;
    SYS.Phi = Phi;
    SYS.omega2 = omega2;
    SYS.plotSettings = plotSettings;
    SYS.modelSettings = modelSettings;
    SYS.simulationSettings = simulationSettings;
    SYS.sys = sys;


%{
The total model looks like:

State:      Size:       Quantity:
q = |z |    Nmodes      Modal states
    |zd|    Nmodes      Modal velocities     

Inputs & Outputs:
Inputs:                           Model:      Outputs:
                               ___________    | L | laser measurement 
| F | Force                   |           |   |ps1| Piezo Sensor outputs
|pa1| Piezo Actuator inputs   |           |   |...|
|...|                       =>|    sys    |=> |psn|        
|pan|                         |           |   |ac1| Accelerometer outputs  
                              |___________|   |...|
                                              |acn|
The amount of piezo actuators and sensors depends on the settings at the
defined top of the file. Also the amount of accelerometers. The system is
by default in modal space (TODO).

SYS contains all necessary information for plotting and
simulation. This way SYS can just be saved and used in other files or
functions. 
%}
   
%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be able to simulate the system, it is discretised. 
dsys = c2d(sys,simulationSettings.dt,'ZOH');            % Discretise using ZOH
dSYS = SYS;                                             % Make new system
dSYS.sys = dsys;
dSYS.descr = 'Discrete time beam model';

%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulates the system and makes nice plots and stuff. It is not
% really necessary and not the point of this script. But worked wonders when
% debugging the code. 

if simulationSettings.simulate ==  true
    disp('Simulating...')

    % Handy defenitions
    T = simulationSettings.T;
    dt = simulationSettings.dt;
    t = 0:dt:T;

    ny = size(dsys,1);                  % Number of outputs
    nu = size(dsys,2);                  % Number of inputs
    nq = size(dsys.A,1);                % Number of states

    % Distubrance input generation (function at the bottom)
    Udist = generateInput(simulationSettings);  % Disturbance input

    % Some initialisations.
    qfull = zeros(nq,length(t));        % Full state vector over time
    y = zeros(ny,length(t));            % System output
    U = zeros(nu,1);                    % Input vector
    q0 = zeros(nq,1);     % Normal time simulation

    % Set up waitbar if necessary
    if simulationSettings.waitBar == true
        f = waitbar(0,'1','Name','Running simulation loop...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(f,'canceling',0);
    end

    % Simulation loop!!!! Here, the system is forward Eulered.
    q1 = q0;    % Initialise at q0 (Looks weird, but is correct. (look at first line of loop)
    for i = 1:T/dt+1
        % Shift one time step
        q = q1; % Previous next state is the current state. (Yes, very deep indeed)

        % Pick input
        U(simulationSettings.distInput) = Udist(i);

        % Propagate discrete time dynamic system
        q1 = dsys.A*q + dsys.B*U;           % Propagate dynamics
        y(:,i) = dsys.C*q + dsys.D*U;       % Measurement equation

        % Also save full states for plotting (in q space, so modal)
        qfull(:,i) = q;

        % Run state estimation algorithms here somewhere

        % Update waitbar and check cancel button
        if simulationSettings.waitBar == true
            waitbar(i/(T/dt+1),f,sprintf('Step %d/%d',i,T/dt+1))
            if getappdata(f,'canceling')
                break
            end
        end
    end

    if simulationSettings.waitBar == true
        delete(f);
    end

    % Save simulation data in the model itself. 
    dSYS.simulationData.qfull = qfull;
    dSYS.simulationData.y = y;
    dSYS.simulationData.t = t;
    dSYS.simulationData.Udist = Udist;

    % Make a nice plot of the simulation!
    if plotSettings.plot == true
       plots = plotter(dSYS);
    end
end

disp('Done!')

%% FUNCTIONS!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get element lengths for adaptive mesh (Taking smart patches into account)
function [Ls,sElements] = getLengths(modelSettings,sBeam)

    nPatches = length(modelSettings.patches);   % Number of patches
    nsElementsP = modelSettings.nsElementsP;    % Number of smart elements per patch
    nsElements = nPatches*nsElementsP;          % Number of smart elements in total
    L = modelSettings.L;

    if ~isempty(modelSettings.patches)
        ngaps = nPatches+1;                         % Number of gaps in theory
        gaps = zeros(1,ngaps);                      
    
        for i = 1:ngaps% Get the length of every gap (between the patches)
            if i == 1                   % The first gap
                gaps(i) = modelSettings.patches(i)*L; 
            elseif i <= nPatches                       % The rest of the gaps
                gaps(i) = modelSettings.patches(i)*L-(sum(gaps(1:i))+(i-1)*sBeam.L*nsElementsP);
            else
                gaps(i) = L-(sum(gaps(1:i))+(i-1)*sBeam.L*nsElementsP);
            end
        end
    
        % Catch some errors here already
        if sum(gaps)+nPatches*sBeam.L*nsElementsP - L > 0
            error('Gaps are not correct (do not sum up to total length of beam)')
        end
        for i = 1:length(gaps)
            if gaps(i) < 0 
                error('Patches are too close together! (They overlap)')
            end
        end
    
        nGapElements = ceil(gaps/(modelSettings.LbElements*L)); % Get number of elements per gap
    
        gelementLengths = cell(ngaps,1);                         % Intermediate cell array containing lengths of elements in gaps
        for i = 1:ngaps
            elementL = gaps(i)/nGapElements(i);
            gelementLengths{i} = ones(1,nGapElements(i))*elementL;
        end
    
        pelementLengths = cell(nPatches,1);
        for i = 1:nPatches
            elementL = sBeam.L;
            pelementLengths{i} = ones(1,modelSettings.nsElementsP)*elementL;
        end
    
        Ls = [];
        sElements = [];
        for i = 1:nPatches+1
           if i <= nPatches
               Ls = [Ls,gelementLengths{i},pelementLengths{i}];
               sElements = [sElements,zeros(1,length(gelementLengths{i})),ones(1,length(pelementLengths{i}))];
           else
               Ls = [Ls,gelementLengths{i}];
               sElements = [sElements,zeros(1,length(gelementLengths{i}))];
           end
        end
        
        if round(sum(Ls) - L,4) > 0
            error('Ls is not correct! (not the same as L')
        end
        if sum(sElements) - nsElements > 0
            error('sElements is incorrect!')
        end
    else
        nElements = 1/modelSettings.LbElements;
        Ls = ones(1,nElements)*L/nElements;
        sElements = [];
    end
end

% Compute analytical integrals for paper. 
function [intS,intG] = shapeFunctions()
    syms L y
    
    A = [1, 0, 0, 0;
            0, 1, 0, 0;
            1, L, L^2, L^3;
            0, 1, 2*L, 3*L^2];
    N = [1, y, y^2, y^3]/A;
    n1 = diff(N,y);
    n2 = diff(n1,y);
    
    intS = int(n2,0,L);
    intG = int(n1,0,L);
end

% Input generation
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
        startSample = ceil(simulationSettings.stepTime(1)/dt);
        uStep(startSample:end) = ones(1,T/dt-startSample+2)*simulationSettings.stepAmp;
        if length(simulationSettings.stepTime) > 1                          % If end time is defined
            stopSample = floor(simulationSettings.stepTime(2)/dt);
            uStep(stopSample:end) = zeros(1,T/dt-stopSample+2);
        end
    end

    % Inpulse input
    uImpulse = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.impulseTime)
        impulseSample = ceil(simulationSettings.impulseTime/dt);
        uImpulse(impulseSample) = simulationSettings.impulseAmp;
    end

    % harmonic input
    uHarmonic = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.harmonicTime)
        startSample = ceil(simulationSettings.harmonicTime/dt);
        t = simulationSettings.harmonicTime-dt:dt:T;
        A = simulationSettings.harmonicAmp;
        f = simulationSettings.harmonicFreq;
        uHarmonic(startSample:end) = A*sin(2*pi*t/f);
    end

    % Random Input
    uRand = zeros(1,T/dt+1);
    if ~isempty(simulationSettings.randTime)
        startSample = ceil(simulationSettings.randTime/dt);
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