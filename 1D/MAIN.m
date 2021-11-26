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
%}

clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

disp('Setting up...')
%% Parameters & Settings
% Basic settings
L = 370e-3;     % Total length
b = 40e-3;      % Total width
patchL = 50e-3; % patch length

% Model settings
modelSettings.patches = []; %[0,0.25];                          % location of start of patches (y/L)
modelSettings.nsElementsP = 3;                          % Number of smart elements per patch
modelSettings.LbElements = 0.1;                         % Preferred length of beam elements (y/L) (will change slightly)

modelSettings.dispInput = false;                        % Use displacement at the base as input?
    modelSettings.kBase = 1e10;                         % Stiffness of spring at base
    modelSettings.cBase = 0.01;
modelSettings.Nmodes = 10;                               % Number of modes to be modelled (Can't be larger then the amount of nodes)
modelSettings.measurementHeight = 1;                 % Height of the measurement

modelSettings.fb = false;                               % Apply feedback control law?

modelSettings.L = L;                                    
modelSettings.b = b;

% plotSettings
plotSettings.plotNodes = true;                          % Plot the nodes
    plotSettings.nodeNumbers = false;                   % Plot the node numbers
plotSettings.elementNumbers = true;                     % Plot the element numbers
plotSettings.sensor = true;                             % Plot the sensor in beam plot (red line)
plotSettings.sensorPlot = true;                         % Plot the sensor measurement

% simulationSettings
simulationSettings.simulate = true;                     % Simulate at all?
    simulationSettings.stepTime = 0.1;                  % Location of input step
    simulationSettings.dt = 1e-3;                       % Time step
    simulationSettings.T = 2;                           % Total simulation time

simulationSettings.plot = true;                         % Plot the simulation?
    simulationSettings.animate = false;                 % Animate the plot?
        simulationSettings.pauseStart = false;          % Pause before animation (To start screen recorder)

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

%% Matrix setup 
nPatches = length(modelSettings.patches);           % Number of patches
nsElementsP = modelSettings.nsElementsP;            % Number of elements

% Get element lengths! (Function at the bottom)
[Ls,sElements] = getLengths(modelSettings,sBeam);   % Vector with for every element its length, and which elements are smart
modelSettings.sElements = find(sElements);          % Which elements are smart?

% Assembly
numEl = length(Ls);
numNodes = numEl+1;
Nmodes = modelSettings.Nmodes;

K = zeros(numNodes*2);
M = zeros(numNodes*2);
elements = element.empty(numEl,0);

% Loops over all elements, constructs their elemental K and M matrices and
% stores them in element objects aswell (class defenition in element.m).
% Also constructs the system K and M matrices by putting elemental matrices
% in there. 
for i = 1:numEl
    % Elemental coordinates: q = [v1,theta1,v2,theta2]';
    n1 = i;     % Starting node
    n2 = i+1;   % Ending node
  
    if sum(ismember(modelSettings.sElements,i)) > 0 % Is smart element?
        El = copy(sBeam);
        El.sBeam = true;
    else                                            % Or normal beam element?
        El = copy(Beam);
    end
    
        El.L = Ls(i);   % Give element correct length
        
        % Calculate stiffness 
        Ib = El.b*El.h^3/12;
        Ip = El.pb*El.ph^3/12 + El.pb*El.ph*(El.ph+El.h)^2/4;
        
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
    
    % Apply boundary condition to first node
    if n1 == 1                                      % If the node we are looking at is the first one
        if modelSettings.dispInput == false         % No displacement input at the base (Only force)
            % Prescribe 0 to displacement and rotation
            Ke(1:2,:) = zeros(2,4);                   
            Ke(:,1:2) = zeros(4,2);
            Ke(1:2,1:2) = eye(2);                   % Fix both rotation and displacement
            
            % Prescribe 0 to acceleration aswell
            Me(1:2,:) = zeros(2,4);
            Me(:,1:2) = zeros(4,2);
            Me(1:2,1:2) = eye(2);
        else                                        
            Ke(2,:) = zeros(1,4);                   
            Ke(:,2) = zeros(4,1);
            Ke(2,2) = 1;                            % Fix only the rotation (not the displacement)

            Me(2,:) = zeros(1,4);                   
            Me(:,2) = zeros(4,1);
            Me(2,2) = 1;                            % Fix only the rotary acceleration (not the displacement)
        end
    end
    
    % Assembly. Since everything is nicely connected in sequence (1to2to3..) , the
    % elemenental matrices stay nicely together.
    K(n1*2-1:n2*2,n1*2-1:n2*2) = K(n1*2-1:n2*2,n1*2-1:n2*2) + Ke;   % K
    M(n1*2-1:n2*2,n1*2-1:n2*2) = M(n1*2-1:n2*2,n1*2-1:n2*2) + Me;   % M
    
    %% Put information in element objects aswell (mainly for plotting, also for structuring)
    x1 = 0;
    x2 = 0;
    y1 = sum(Ls(1:i-1));
    y2 = sum(Ls(1:i));
    t1 = 0;
    t2 = 0;
    initPos = [x1, x2;
                y1, y2;
                t1, t2];
    
    elements(i) = El.assign(i,[n1,n2],initPos,Ke,Me);
end

% Find mode shapes and natural frequencies
if modelSettings.Nmodes <= numEl
    if modelSettings.dispInput == true
%         [Phi,omega2] = eigs(K,M,modelSettings.Nmodes+1,'smallestabs');      % Compensate for 1 weird mode
        [Phi,omega2] = eig(K,M);
    else
        [Phi,omega2] = eig(K,M);
    end
else
    error('Number of elements is smaller then the number of modes requested. Modal decomposition will not work!')
end

Phi = Phi(:,1:Nmodes+2);
omega2 = omega2(1:Nmodes+2,1:Nmodes+2);

%Snip off weird modes 
%{
    When displacement input is turned on, the first eignemode is a rotation
    around the first node. This should not be possible as this rotation is
    fixed. The other modes are correct however. 

    When dispInput is turned off, the first node is fixed but the first two
    eigenmodes do show a movement of this node. This should not be
    possible. The other modes are correct however. 
    
    I think this is not very influential, they are snipped off, but might
    need to look at this later if model turns out inaccurate. 
%}
if modelSettings.dispInput == true
    Phi = Phi(:,3:end);
    omega2 = omega2(3:end,3:end);
else 
    Phi = Phi(:,3:end);
    omega2 = omega2(3:end,3:end);  
end

% Check eigenfrequencies
f = diag(sqrt(omega2)/(2*pi));

I = (Beam.h^3)/12;
S = Beam.h;         % Need to check this formulation still!
lambda = 1.87510407;
analyticalf1 = lambda^2/(2*pi*L^2)*sqrt(Beam.E*I/(Beam.rho*S));

delta = 1*L^3/(3*Beam.E*I);

if abs(f(1) - analyticalf1) > 0.001 && nPatches == 0 && modelSettings.dispInput == false  % Check only when there are no patches
    warning('Eigenfrequencies do not correspond')
end

%% Dynamical matrices setup
Phi = Phi/(Phi'*M*Phi);     % Normalise w.r.t. mass matrix

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
% patches
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

% Damping matrix 
Cmodal = 2*Beam.zeta*sqrt(omega2);

% External input B matrix (in d coordinates)
Bext = zeros(numNodes*2,1);
Bext(end-1) = 1;                          % Force at last node in x direction 

% Base reference B matrix (in d coordinates)
Bx = zeros(numNodes*2,1);
if modelSettings.dispInput == true
    Bx(1) = modelSettings.kBase;
end

Bxd = zeros(numNodes*2,1);
if modelSettings.dispInput == true
    Bxd(1) = modelSettings.cBase;
end

% Base location C marix (in d coordinates)
Cx = zeros(1,numNodes*2);
Cx(1) = 1;

% Measurement C matrix (in d coordinates) (have to interpolate between two
% nodes)
    % Find interpolate nodes
height = modelSettings.measurementHeight*L;
pos = unique([elements.pos]');
diff = pos-height;
lowerNodes = find(diff<0);
interpNodes = [lowerNodes(end),lowerNodes(end)+1];  

    % Corresponding elements
for i = 1:numEl
    if interpNodes == elements(i).neighbours          
        interpEl = i;
    elseif interpNodes(1) == elements(i).neighbours(2) && interpNodes(2) > numNodes % If exactly at last node
        interpEl = i-1;
        interpNodes = interpNodes-1;
    end
end

    % natural coordinate alpha
alpha = height - pos(interpNodes(1));
N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
Cmeasurement = zeros(1,numNodes*2);                         
Cmeasurement(interpNodes(1)*2-1:interpNodes(1)*2+2) = N;  % Construct C from this (cubic interpolation)

% Construct full system 
A = [zeros(Nmodes),eye(Nmodes);
    -omega2,-Cmodal];
B = [zeros(Nmodes,nPatches+3);                          
    Phi'*Bext, Phi'*Bx, Phi'*Bxd, Phi'*Bg];                       % [Force input, Base force, Base force, Piezo inputs]
C = [Cmeasurement*Phi, zeros(1,Nmodes);                 % Laser measurement
    Cx*Phi,zeros(1,Nmodes);                             % Base position
    zeros(1,Nmodes),Cx*Phi;                             % Base velocity
    Cs*Phi,zeros(nPatches,Nmodes)];                     % Piezo outputs

sys = ss(A,B,C,[]);
if modelSettings.dispInput == true          
    sys = feedback(sys,eye(2),[2,3],[2,3],-1);
end

% Cut out unescessary inputs and outputs!
tvec = 1:3+nPatches;
tvec(3) = [];
sys = sys(tvec,tvec);

%{
The total model looks like:

State:      Size:       Quantity:
q = |z |    Nmodes      Modal states
    |zd|    Nmodes      Modal velocities     

Inputs & Outputs:
Inputs          Model           Outputs
|F |Force                       |L | laser measurement 
|r |Base Ref     ___________    |d | base Location 
|a1|Act inputs  |           |   |s1| Sensor outputs
|a2|            |           |   |s2|
|a3|          =>|   sys     |=> |s3|        
|..|            |           |   |..|   
|an|            |___________|   |sn|
%}

% Apply feedback
% if there are no patches, no feedback can be done
if nsElements == 0
    modelSettings.fb = false;    
end

if modelSettings.fb == true
    [~,wpeak] = hinfnorm(sys(1,1));
    wc = wpeak;
    zetac = 10*Beam.zeta;
    kc = 100;

    PPF = tf(kc*wc^2,[1 2*zetac*wc wc^2]);
    PPF = PPF*eye(nPatches);                                    % SISO PPF for every patch (just to test)

    sys_fb = feedback(sys,PPF,[3:2+nPatches],[3:2+nPatches],1); % Same system as before, only with ppf
end    

%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulates the system and makes nice plots and stuff. It is not
% really necessary and not the point of this script. But worked wonders whe
%n debugging the code. 

if simulationSettings.simulate ==  true
    disp('Simulating...')
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    
    q0sim = zeros(size(sys.A,1),1);     % Normal time simulation
    Asim = sys.A;
    Bsim = sys.B;
    [t,q] = ode23t(@(t,z) sysFun(t,z,Asim,Bsim,modelSettings,simulationSettings),tvec,q0sim);
    q = q';
    y = sys.C*q;
    dmat = Phi*q(1:Nmodes,:);
    measurement = y(1,:);
    baseLocation = y(2,:);
    
    if modelSettings.fb == true
        disp('Simulating fb...')
        q0sim = zeros(size(sys_fb.A,1),1);  % Feedback time simulation
        Afb = sys_fb.A;
        Bfb = sys_fb.B;
        [t_fb,q_fb] = ode45(@(t,z) sysFun(t,z,Afb,Bfb,modelSettings,simulationSettings),tvec,q0sim);
        q_fb = q_fb';
        y_fb = sys_fb.C*q_fb;
        measurement_fb = y_fb(1,:);
    end
    
    if simulationSettings.plot == true
        disp('Plotting simulation...')
        figure()
        %sgtitle(['Beam with ',modelSettings.Input,' input'])
        subplot(7,3,[1,4,7,10]) % Beam plot
            hold on
            grid on
            xlabel 'm'
            ylabel 'm'
            axis equal
            xlim([-L/6,L/6]);
            ylim([0,1.2*L])
            title 'Beam plot'
            beamAx = gca;
        subplot(7,3,[2,3,5,6])
            hold on
            grid on
            ylabel 'displasement [m]'
            title 'Laser Measurement'
            measurementAx = gca;
            xlim([0,simulationSettings.T]);
        subplot(7,3,[8,9,11,12]);
            hold on
            grid on
            xlabel 'time [s]'
            ylabel 'V'
            title 'Output voltages'   
            stateAx = gca;
            xlim([0,simulationSettings.T]);
        subplot(7,3,[13:21])
            bodeAx = gca;
            
        % Plot bode
        plotoptions = bodeoptions;
        plotoptions.Title.String = '';
        plotoptions.Title.Interpreter = 'latex';
        plotoptions.XLabel.Interpreter = 'latex';
        plotoptions.YLabel.Interpreter = 'latex';
        plotoptions.XLabel.FontSize = 9;
        plotoptions.YLabel.FontSize = 9;
        plotoptions.FreqUnits = 'Hz';
        plotOptions.grid = 'on';
        plotOptions.PhaseWrapping = 'off';
        plotOptions.XLimMode = {'manual','manual'};
  
        if modelSettings.fb == true
            bod = bodeplot(bodeAx,sys(1,1),sys_fb(1,1),plotoptions);
        else
            bod = bodeplot(bodeAx,sys(1,1),plotoptions);
        end
        
        % Plot analytical omega1
        xline(bodeAx,analyticalf1)
      
        simPlots = [];     
        
        if simulationSettings.animate == true
            Nsteps = length(t);
        else
            Nsteps = 1;
        end

        % Run animation loop! (Only once if animate is turned off)
        for i = 1:Nsteps
            tstart = tic;
            if simulationSettings.animate == true
                d = dmat(:,i);
            else
                d = dmat(:,end);
            end

            % update elements
            for j = 1:length(elements)
                elements(j).update(d);
            end
            
            % Delete old plot
            if ~isempty(simPlots)
                delete(simPlots)
            end

            % Plot beam
            for j = 1:numEl
                pel = elements(j).show(beamAx,plotSettings); 
                simPlots = [simPlots,pel];
            end
            
            timeText = text(beamAx, 0,L*1.1,['Time: ',num2str(round(t(i),1)),'/',num2str(t(end)),' s'],...
                                                        'HorizontalAlignment','center');
            simPlots = [simPlots, timeText];
            
            % Plot laser
            if plotSettings.sensor == true
                if simulationSettings.animate == true
                    laserx = measurement(i);
                else
                    laserx = measurement(end);
                end
                laser = plot(beamAx,[beamAx.XLim(1) laserx],[1,1]*modelSettings.measurementHeight*L,'r','lineWidth',2);
                simPlots = [simPlots, laser];
            end
            
            % Plot laser measurement & Voltages
            if plotSettings.sensorPlot == true
                if simulationSettings.animate == true
                    if i > 1
                        xs = [t(i-1),t(i)];
                        ys = [measurement(i-1),measurement(i)];
                        if nsElements > 0               % If there are patches
                            ysbeam = [y(3:nPatches+2,i-1),y(2:nPatches+1,i)];
                            if modelSettings.fb == true
                                y_fb = [measurement_fb(i-1),measurement_fb(i)];
                            end
                        end
                    else
                        xs = t(i);
                        ys = measurement(i);
                        if nsElements > 0               % If there are patches
                            ysbeam = y(3:nPatches+2,i);
                            if modelSettings.fb == true
                                ys_fb = measurement_fb(i);
                            end
                        end
                    end
                else
                    xs = t;
                    baseLocation = dmat(1,:);
                    ys = measurement;
                    if nsElements > 0                   % If there are patches
                        ysbeam = y(3:nPatches+2,:);
                        if modelSettings.fb == true      
                            ys_fb = measurement_fb;
                        end
                    end
                end
                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
                    if nsElements > 0                   % If there are patches
                        statePlot = plot(stateAx,xs,ysbeam);
                        if modelSettings.fb == true
                            sensorfbPlot = plot(measurementAx,xs,ys_fb,'color',[0.8500 0.3250 0.0980]);
                        end
                    end
            end
            
            drawnow;
            
            % Measure elapsed time and delay
            elapsed = toc(tstart);
            
            if simulationSettings.animate == true
                %pause(max(simulationSettings.dt-elapsed,0));
            end
            if simulationSettings.pauseStart == true && i == 1 && simulationSettings.animate == true
                disp('PAUSED, Press any key to continue');
                pause;
            end
        end 
    end
end

disp('Done!')


% FUNCTIONS!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMode(n,Phi,elements,plotSettings)
    numEl = length(elements);
    d = Phi(:,n);
    d = d/max(abs(d));
    
    localElements = copy(elements);
    
    figure()
    hold on
    axis equal
    for i = 1:numEl
        localElements(i).update(d);
        localElements(i).show([],plotSettings);
    end
end

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
    numEl = sum(nGapElements)+nsElements;                   % Total number of elements

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

% odefun for ode45 etc..
function qd = sysFun(t,q,A,B,modelSettings,simulationSettings)
    % Generic input!
    nsElements = length(modelSettings.sElements);
    if t >= simulationSettings.stepTime
        F = 1;
        r = 0;%5e-3;
        V = 0;    % Input voltage on the patch (I think)
    else
        F = 0;
        r = 0;
        V = 0;
    end
    nPatches = length(modelSettings.patches);
    N = ones(1,nPatches);
    % N(4) = 1; % Which patch?
    % Evaluate system
    U = [F,r,N*V]';
    qd = A*q+B*U;

    if sum(isnan(q))>0 || sum(isnan(qd)) > 0
        error("NaN's found in state vector")
    end        
end 