function [SYS,modelSettings] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS)
%% Matrix setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPatches = length(modelSettings.patches);           % Number of patches
nsElementsP = modelSettings.nsElementsP;            % Number of elements per patch
nAcc = length(modelSettings.Acc);                   % Number of accelerometers
L = modelSettings.L;

% Get element lengths! (Function at the bottom)
[Ls,sElements] = getLengths(modelSettings,sBeam);   % Vector with for every element its length, and which elements are smart
modelSettings.sElements = sElements;          % Which elements are smart?

% Find accelerometer locations
accHeights = L*modelSettings.Acc;
nodeHeights = [0,cumsum(Ls)];
accPos = zeros(nAcc,2);         % accelerometer positions [Element, eta]

for i = 1:nAcc
    accHeight = accHeights(i);
    difference = nodeHeights-accHeight;
    upperNode = find(round(difference,4)>=0,1);           % below this node
    lowerNode = upperNode-1;                        % Above this node
    accElement = lowerNode;                         % This element!
    
    if lowerNode == 0 % If accelerometer is placed at 0. 
        lowerNode = 1;
        upperNode = 2;
        accElement = 1;
        eta = 0;
    else
        % Natural element coordinate (How far along this element, for integration)
        eta = (accHeight-nodeHeights(lowerNode))/(nodeHeights(upperNode)-nodeHeights(lowerNode));
    end

    accPos(i,1) = accElement;
    accPos(i,2)= eta;
    accPos(i,3)= i;
end

modelSettings.accElements = accPos;                 % Which elements have Accelerometers?

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
  
    if any(ismember(modelSettings.sElements,i),'all')     % Is smart element?
        El = copy(sBeam);                           % Inherit all smart element parameters
        El.sBeam = true;                            % Tell it it's a smart element
        El.piezoNumber = find(any(ismember(modelSettings.sElements,i)));
        El.piezoElements = modelSettings.sElements(:,El.piezoNumber)';
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
    if n1 == 1
        % Prescribe 0 to displacement and rotation
        Ke(1:2,:) = zeros(2,4);                   
        Ke(:,1:2) = zeros(4,2);
        Ke(1:2,1:2) = eye(2);                   % Fix both rotation and displacement
        
        % Prescribe 0 to acceleration aswell
        Me(1:2,:) = zeros(2,4);
        Me(:,1:2) = zeros(4,2);
        Me(1:2,1:2) = eye(2);
    end

    % Apply boundary condition to last node
    if n2 == numNodes                                      % If the node we are looking at is the first one
        % Prescribe 0 to displacement and rotation
        Ke(4,:) = zeros(1,4);                   
        Ke(:,4) = zeros(4,1);
        Ke(4,4) = 1;                   % Fix only rotation
        
        % Prescribe 0 to acceleration aswell
        Me(4,:) = zeros(1,4);
        Me(:,4) = zeros(4,1);
        Me(4,4) = 1;
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

Phi = Phi(:,1:Nmodes+3);                % Select only Nmodes+2 modes
omega2 = omega2(1:Nmodes+3,1:Nmodes+3); % Select only Nmodes+2 natural frequencies

%Snip off weird modes 
Phi = Phi(:,4:end);                     % First 2 modes are not correct, so snip them off
omega2 = omega2(4:end,4:end);           % Same for natural frequencies

% Check first natural frequency
f = diag(sqrt(omega2)/(2*pi));
I = (Beam.h^3)/12;
S = Beam.h;         
lambda = 1.87510407;
analyticalf1 = lambda^2/(2*pi*L^2)*sqrt(Beam.E*I/(Beam.rho*S));
  
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
        height = modelSettings.forceHeight*L;
        pos = unique([elements.pos]');
        difference = pos-height;
        upperNode = find(round(difference,4)>=0,1);
        interpNodes = [upperNode-1,upperNode];  
        interpEl = upperNode-1;

        % local coordinate alpha
        alpha = height - pos(interpNodes(1));                       % alpha = eta*obj.L
        N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
        
        Bext = zeros(numNodes*2,1);
        Bext(interpNodes(1)*2-1:interpNodes(1)*2+2) = N;
        
% C matrices    
    % Measurement C matrix for laser measurement
        % Find interpolate nodes
        height = modelSettings.measurementHeight*L;
        pos = unique([elements.pos]');
        difference = pos-height;
        upperNode = find(round(difference,4)>=0,1);
        interpNodes = [upperNode-1,upperNode];  
        interpEl = upperNode-1;
        
        % local coordinate alpha
        alpha = height - pos(interpNodes(1));                       % alpha = eta*obj.L
        N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
        Cmeasurement = zeros(1,numNodes*2);                         
        Cmeasurement(interpNodes(1)*2-1:interpNodes(1)*2+2) = N;    % Construct C from this (cubic interpolation)
    
    % Acceleration measurement C matrix
        % Accelerometer C matrix
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

if modelSettings.strainRate == true
    C = [Cmeasurement*Phi, zeros(1,Nmodes);                 % Laser measurement
        zeros(nPatches,Nmodes),Cs*Phi;                      % Piezo outputs
        -Cacc*Phi*omega2, -Cacc*Phi*Cmodal];                % Accelerometer outputs
else
    C = [Cmeasurement*Phi, zeros(1,Nmodes);                 % Laser measurement
    Cs*Phi,zeros(nPatches,Nmodes);                      % Piezo outputs
    -Cacc*Phi*omega2, -Cacc*Phi*Cmodal];                % Accelerometer outputs
end

D = [zeros(nPatches+1,size(B,2));                       % No throughput laser,base and piezo measurements
    Cacc*(Phi*Phi')*[Bext,Bg]];                  % Throughput for acceleration measurements!

sys = ss(A,B,C,D);

% Define true process and measurement noise covariances
Q = eye(Nmodes*2)*modelSettings.wcov;
% R = [modelSettings.laserCov,zeros(1,nPatches+nAcc);
%     zeros(nPatches,1), eye(nPatches)*modelSettings.patchCov, zeros(nPatches,nAcc);            % Allow for different covariances for different sensors!
%     zeros(nAcc,1), zeros(nAcc,nPatches), eye(nAcc)*modelSettings.accCov]; 
R = blkdiag(modelSettings.laserCov, eye(nPatches).*modelSettings.patchCov, eye(nAcc)*modelSettings.accCov);
S = zeros(Nmodes*2,1+nPatches+nAcc);


% Noise influence matrices Bw and Dv
Bw = eye(Nmodes*2);
Dv = eye(1+nPatches+nAcc);

% Put the system together with other information in a larger model object. 
SYS.elements = elements;
SYS.analyticalf1 = analyticalf1;
SYS.numEl = numEl;
SYS.numNodes = numNodes;
SYS.Nmodes = Nmodes;
SYS.Phi = Phi;
SYS.omega2 = omega2;
SYS.plotSettings = plotSettings;
SYS.modelSettings = modelSettings;
SYS.simulationSettings = simulationSettings;
SYS.Q = Q;
SYS.R = R;
SYS.S = S;
SYS.Bw = Bw;
SYS.Dv = Dv;
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
                              |           |   |...|
                            =>|    sys    |=> |psn|        
                              |           |   |ac1| Accelerometer outputs  
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
dsys = c2d(sys,simulationSettings.dt,modelSettings.c2dMethod);            % Discretise using ZOH
SYS.dsys_sim = dsys;                                                 % Make new system

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    sElements = find(sElements);
    sElements = reshape(sElements,modelSettings.nsElementsP,[]);
end

% Compute analytical integrals for paper. 
function [intS,intG] = shapeFunctions()
    y = [];
    Lss = [];
    syms Lss y
    
    Ainterp = [1, 0, 0, 0;
            0, 1, 0, 0;
            1, Lss, Lss^2, Lss^3;
            0, 1, 2*Lss, 3*Lss^2];
    shapeFunction = [1, y, y^2, y^3]/Ainterp;
    dshapeFunction = diff(shapeFunction,y);
    ddshapeFunction = diff(shapeFunction,y);
    
    intS = int(ddshapeFunction,0,Lss);
    intG = int(dshapeFunction,0,Lss);
end

end