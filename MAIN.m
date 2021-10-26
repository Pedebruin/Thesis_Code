% This is the main file.

clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

% Create template objects
beam = body('Beam');
sensor = body('Sensor');
actuator = body('Actuator');

%% Parameters & Settings
% beam parameters
beam.L = 370e-3;                    % m
beam.h = 2e-3;                      % m
beam.b = 40e-3;                     % m
beam.N = 65;%L/h;                   % Number of Vertical Nodes (minimum 2) (Accurate from around 65)
beam.E = 70E9;                      % E modulus
beam.G = 25.5E9;                    % Shear Modulus
beam.mu = 0.334;                    % Poisson
beam.rho = 2710;                    % Mass density
beam.alphaC = 0.001599;             % For proportional damping
beam.betaC = 0.001595;              % For proportional damping
beam.zeta = 0.01;                   % For modal damping

% sensor parameters
simulationSettings.Nsensors = 5;    % number of sensor patches
sensor.L = 30e-3;                   % m
sensor.h = 1e-3; %1;                % m
sensor.b = 40e-3;                   % m
sensor.N = 3;
sensor.E = 5.2e10;
sensor.G = 0.775e9;                 % Shear modulus
sensor.mu = 0.334;
sensor.rho = 7800;
sensor.alphaC = 0.001599;
sensor.betaC = 0.001595;
sensor.zeta = 0.01;

sensor.d31 = 2.2e-11;
sensor.s11 = -3.0e-11;
sensor.eta33 = 2400*sensor.eta0;

% actuator parameters
simulationSettings.Nactuators = 3;  % number of actuators
actuator.L = 30e-3;                 % m
actuator.h = 1e-3; %1               % m
actuator.b = 40e-3;                 % m
actuator.N = 3;
actuator.E = 5.2e10;
actuator.G = 0.775e9;               % Shear modulus
actuator.mu = 0.334;
actuator.rho = 7800;
actuator.alphaC = 0.001599;
actuator.betaC = 0.001595;
actuator.zeta = 0.01;

actuator.d31 = 2.2e-1;
actuator.s11 = -3.0e-1;
actuator.eta33 = 2400*actuator.eta0;

wmax = 1e4*(2*pi);                  % Maximum frequency to include in modal decomposition

% Simulation settings
simulationSettings.simulate = true;
    simulationSettings.plot = true;
    simulationSettings.animate = false;
        simulationSettings.pauseStart = false;
    simulationSettings.stepTime = 0.5;              % s
    simulationSettings.T = 5;                       % s
    simulationSettings.dt = 0.01;                   % s 

% model settings
simulationSettings.measurementHeight = 0.85;    % []*L Height of laser measurement
simulationSettings.forceHeight = 1;             % []*L Height of applied force
simulationSettings.dampingModel = 'modal';      % 'Modal','proportional' or 'none'
simulationSettings.Input = 'disp';             % 'force' for force input, or 'disp' for disp input   
    simulationSettings.Kbase = 1e10;            % only in 'disp' mode

% Plot settings
% Single mode plot settings
plotSettings.modalAnalysis = false;          % Single plot with a specific mode
    plotSettings.modeNr = 2;
    plotSettings.A = 0.5;
    
% Beam, measurement and bode plot
plotSettings.realMesh = false;                  % Make an extra plot with the original Mesh
plotSettings.color = 'none';                    % 'disp' for displacement or 'deform' for deformation or 'none'
    plotSettings.fill = true;                   % Fill the beam with patch
plotSettings.nodes = false;                     % Plot the nodes 
    plotSettings.nodeNumbers = false;           % Node numbers
plotSettings.links = true;                     % Plot the links (edges)
    plotSettings.linkNumbers = false;           % link numbers
plotSettings.sensor = true;                     % plot the sensor
plotSettings.sensorPlot = true;                 % plot the sensor output
plotSettings.base = true;                       % Plot the location of the base in the measurement


%% Build FEM model of the beam
disp('Running!-------------------------------------------------------------')
[beam,sensors,actuators,FEM,PDEbeam] = buildModel(beam,sensor,actuator,plotSettings,simulationSettings);
bodies = [beam,sensors,actuators];

numBodies = length(bodies);
numNodes = length(unique([bodies.nodes]));
numElements = length(unique([bodies.elements]));
numLinks = length(unique([bodies.links]));

K = FEM.K;
M = FEM.M;

%% Setup LTI model
% Find node to sense!
measureNodes = interpolate(bodies,simulationSettings.measurementHeight*beam.L);
simulationSettings.measurePoint = measureNodes;
                                    
% Find nodes for force input
forceNodes = interpolate(bodies,simulationSettings.forceHeight*beam.L);
simulationSettings.forcePoint = forceNodes;

% Find nodes for displacement input
beamNodes = [beam.nodes.pos];
baseNodes = find(beamNodes(2,:) == 0);
baseNodes = [beam.nodes(baseNodes).number]';
                                                                
% Setup Model
% Input
Bd = zeros(numNodes*2,2);                                       % External Input matrix in d space

forceAlpha = forceNodes(1,4);
Bd(forceNodes(:,1),1) = [1-forceAlpha;forceAlpha];              % Force input (in x)!
Bd(baseNodes,2) = 1/length(baseNodes);                          % For displacement input later on (in x)

% Output
Cd = zeros(2,numNodes*2);                                       % Measurement matrix in d space

measureAlpha = measureNodes(1,4);
Cd(1,measureNodes(:,1)) = [1-measureAlpha;measureAlpha];                      % Measure at selected point
Cd(2,baseNodes) = 1/length(baseNodes);                          % Measure base displacement

% modal system modelling
% Get modal results
wmin = -0.1;    % -0.1 rad/s to capture everything from 0 ->

% Modal decomposition! P395
modal = solve(PDEbeam,'FrequencyRange',[wmin,wmax]);   % Solve
Nmodes = length(modal.NaturalFrequencies);

Phi = [modal.ModeShapes.ux;                                     % Phi
    modal.ModeShapes.uy];
Phi = Phi/(Phi'*M*Phi);

omegaSq = Phi'*K*Phi;                                           % Omega^2

% Damping matrix C_Phi
switch simulationSettings.dampingModel                          % Damping matic C_phi
    case {'modal'}
        C_phi = diag(2*beam.zeta*modal.NaturalFrequencies);
    case {'proportional'}
        C_phi = beam.alphaC*eye(Nmodes)+beam.betaC*omegaSq;
    case {'none'}
        C_phi = zeros(Nmodes);
end
 
Rphi = Phi'*Bd;                                                     % Modal input matrix
Cphi = Cd*Phi;                                                      % Modal measurement

% Assemble new system for modal
disp('Assembling LTI model...')
A_modal = [zeros(Nmodes), eye(Nmodes);                              % Dynamics
     -omegaSq,-C_phi];
B_modal = [zeros(Nmodes,2);Rphi];                                   % External input
C_modal = [Cphi,zeros(2,Nmodes)];                                   % Measurement
D_modal = 0;                                                        % Assume no direct feedthrough
SS_modal = ss(A_modal,B_modal,C_modal,D_modal);                     % 2x2 system!
SS_modal.OutputName = {'Sensor out','Base out'};
SS_modal.Inputname = {'Force in','Base in'};


% Apply stiff spring to base  
switch simulationSettings.Input
    case {'force'}
        % The base is fixed, so no input can be given there. The input is
        % thus chosen as a force at the tip. in x direction.
        SS_modalr = SS_modal(1,[1,2]);
        
    case {'disp'}
        disp('Adding spring to base')
        % To get the input of the system bo be the displacement of the base
        % of the beam, it needs to be fixed by a stiff spring. This is done
        % by applying state feedback with a gain representing a spring
        % constant. 

        K_spring = ss(simulationSettings.Kbase);    % get stiffness matrix

        SS_modal = feedback(SS_modal,K_spring,2,2,-1);
        DCgain = dcgain(SS_modal(1,2));
        SS_modal.B = [B_modal(:,1),B_modal(:,2)/DCgain];

        SS_modalr = SS_modal(1,[1,2]);               % From [F,r] to sensor out. 
        SS_modalr.Inputname = {'Force in','base Ref'};
end

% Check first eigenfrequency!
disp('Checking model...')
I = (beam.h^3)/12;
A = beam.h;
lambda = 1.87510407;
analyticalOmega1 = lambda^2/(2*pi*beam.L^2)*sqrt(beam.E*I/(beam.rho*A));
switch simulationSettings.Input
    case {'force'}
        errPercentage = abs(analyticalOmega1 - modal.NaturalFrequencies(1)/(2*pi))/analyticalOmega1*100;
        
    case {'disp'}
        [~,fpeak] = hinfnorm(SS_modal(1,1));
        errPercentage = abs(analyticalOmega1 - fpeak/(2*pi))/analyticalOmega1*100;
end
if errPercentage > 0.1 % greater then [] procent error in first eigenfrequency!
    warning(['First eigenfrequency does not correspond, use more nodes!',newline...
            'err = ',num2str(round(errPercentage,3)),'%!']);
end

% Augment with Piezo electrics!
disp('Augmenting with piezo states...');
Nsensors = simulationSettings.Nsensors;
Nactuators = simulationSettings.Nactuators;
Np = Nsensors+Nactuators;

if Np > 0
    % Mechanical system matrices q = [z;zd]; u = [F;r], y = [m]
    Amech = SS_modalr.A;
    Bmech = SS_modalr.B;
    Cmech = SS_modalr.C;
    Dmech = SS_modalr.D;

    % For piezo augmentation, the states need to be augmented-> q = [z;zd;q;qd]
    % Where q = [q1,q2...qn] is the charge vector for all piezo patches. 
    patches = [sensors, actuators];


    % Piezo parameters
    [e33s,eta33s] = sensor.setupPiezo(bodies);
    [e33a,eta33a] = actuator.setupPiezo(bodies);
    ds = [0, 0, 0;       % Piezo coupling matrix for PZT in 2D (NEED TO CHECK!)(Depends on direction of polarisation)
        0, 0, 23.3]*10^5;
    etas = eye(2)*eta33s;

    da = [0, 0, 0;       % Piezo coupling matrix for PZT in 2D (NEED TO CHECK!)(Depends on direction of polarisation)
        0, 0, 23.3]*10^5;
    etaa = eye(2)*eta33a;

    T_s = sparse(zeros(Nsensors+Nactuators,numNodes));

    for i = 1:Np
        elements = [patches(i).elements];
        Kpu_s = sparse(zeros(numNodes,numNodes*2));
        Kpp_s = sparse(zeros(numNodes,numNodes));

        switch patches(i).name([1:3])
            case {'Sen'}
                d = ds;
                eta = etas;
            case {'Act'}
                d = da;
                eta = etaa;
        end

        for j = 1:length(elements)
                % Elemental Kpu and Kpp matrices (Nguyen, 2018)
            Kpu_e = sparse(zeros(numNodes,numNodes*2));
            Kpp_e = sparse(zeros(numNodes,numNodes));

                % Vertices
            n1 = elements(j).neighbours(1); % Assume linear elements (Quadratic is slightly more envolved)
            n2 = elements(j).neighbours(2);
            n3 = elements(j).neighbours(3);

                % Positions of the vertices
            x1 = elements(j).pos(1,1);
            x2 = elements(j).pos(1,2);
            x3 = elements(j).pos(1,3);
            y1 = elements(j).pos(2,1);
            y2 = elements(j).pos(2,2);
            y3 = elements(j).pos(2,3);

            A = 1/2*det([1, x1, y1;         % Area of the triangle
                         1, x2, y2;
                         1, x3, y3]);
            Bu = 1/(2*A)*[  y2-y3,  0,      y3-y1,  0,      y1-y2,  0;
                            0,      x3-x2,  0,      x1-x3,  0,      x2-x1;
                            x3-x2,  y2-y3,  x1-x3,  y1-y3,  x2-x1,  y1-y2];
            Bphi = -1/(2*A)*[  y2-y3,  y3-y1,  y1-y2;
                               x3-x2,  x1-x3,  x2-x1];

            Kpu = A*Bphi'*d*Bu;   % Need to look at this ds to derive for 2D!! ('omitted now)
            Kpp = A*Bphi'*eta*Bphi;

            %% Put thes matrices at the right location in the larger elemental
            % matrix!
            Kpux = Kpu(:,[1,3,5]);
            Kpuy = Kpu(:,[2,4,6]);
            Kpu = [Kpux, Kpuy];

            % Kpu matrix
            Kpu_e(n1,n1) = Kpux(1,1);   Kpu_e(n1,n1+numNodes) = Kpuy(1,1);
            Kpu_e(n1,n2) = Kpux(1,2);   Kpu_e(n1,n2+numNodes) = Kpuy(1,2);
            Kpu_e(n1,n3) = Kpux(1,3);   Kpu_e(n1,n3+numNodes) = Kpuy(1,3);
            Kpu_e(n2,n1) = Kpux(2,1);   Kpu_e(n2,n1+numNodes) = Kpuy(2,1);
            Kpu_e(n2,n2) = Kpux(2,2);   Kpu_e(n2,n2+numNodes) = Kpuy(2,2);
            Kpu_e(n2,n3) = Kpux(2,3);   Kpu_e(n2,n3+numNodes) = Kpuy(2,3);
            Kpu_e(n3,n1) = Kpux(3,1);   Kpu_e(n3,n1+numNodes) = Kpuy(3,1); 
            Kpu_e(n3,n2) = Kpux(3,2);   Kpu_e(n3,n2+numNodes) = Kpuy(3,2);
            Kpu_e(n3,n3) = Kpux(3,3);   Kpu_e(n3,n3+numNodes) = Kpuy(3,3); 

            % Kpp matrix
            Kpp_e(n1,n1) = Kpp(1,1);
            Kpp_e(n1,n2) = Kpp(1,2);
            Kpp_e(n1,n3) = Kpp(1,3);
            Kpp_e(n2,n1) = Kpp(2,1);
            Kpp_e(n2,n2) = Kpp(2,2);
            Kpp_e(n2,n3) = Kpp(2,3);
            Kpp_e(n3,n1) = Kpp(3,1);
            Kpp_e(n3,n2) = Kpp(3,2);
            Kpp_e(n3,n3) = Kpp(3,3);

            % Put matrix in larger matrix!
            Kpu_s = Kpu_s+Kpu_e;
            Kpp_s = Kpp_s+Kpp_e;
            T_s(i,[n1,n2,n3]) = [1,1,1];
        end
    end

    % Add all potentials on the nodes for the actuators and the sensors to
    % obtain the potential state of each patch! (Thomas, 2012)
    Kc = Kpu_s';
    Ke = Kpp_s;
    T_s = T_s';

    % actuator input and output vectors
    Bact = [zeros(Nsensors,Nactuators);
            eye(Nactuators)];
    Csens = [eye(Nsensors),zeros(Nsensors,Nactuators),zeros(Nsensors,Np)];

    A = [Amech,                                 [zeros(Nmodes,Np),  zeros(Nmodes,Np);
                                                 -Phi'*Kc*T_s,          zeros(Nmodes,Np)];   
        zeros(Np,Nmodes),   zeros(Np,Nmodes),   zeros(Np,Np),       eye(Np);
        -T_s'*Kc'*Phi,            zeros(Np,Nmodes),   -T_s'*Ke*T_s,                 zeros(Np)];
    B = [Bmech,                     zeros(size(Bmech,1),Nactuators);
        zeros(Np,size(Bmech,2)),    zeros(Np,Nactuators);
        zeros(Np,size(Bmech,2)),    Bact]; 
    C = [Cmech,                     zeros(size(Cmech,1),Np*2);
         zeros(Nsensors,Nmodes*2),  Csens];
    D = zeros(size(Cmech,1)+Nsensors,size(Bmech,2)+Nactuators);

    SS_tot = ss(full(A),full(B),C,D);
    
else
    SS_tot = SS_modalr;
    
end

%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if simulationSettings.simulate ==  true
    disp('Simulating...')
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    
    % Time simulation!
    q0 = zeros(Nmodes*2+Np*2,1);
    Asim = SS_tot.A;
    Bsim = SS_tot.B;
    [t,y] = ode45(@(t,z) sysFun(t,z,Asim,Bsim,simulationSettings),tvec,q0);
    z = y';
    output = SS_tot.C*z;
    
    Cbase = SS_modal.C*z(1:Nmodes*2,:);
    baseLocation = Cbase(2,:);
    measurement = output(1,:);


   switch simulationSettings.Input
        case {'force'}
            plotSettings.base = false;
        case {'disp'}
    end
    
    if simulationSettings.plot == true
        disp('Plotting simulation...')
        figure()
        %sgtitle(['Beam with ',simulationSettings.Input,' input'])
        subplot(4,3,[1,4,7,10])
            hold on
            grid on
            xlabel 'm'
            ylabel 'm'
            axis equal
            xlim([-beam.L/6,beam.L/6]);
%             xlim([-beam.h*4,beam.h*4])
            ylim([0,1.2*beam.L])
            title 'Beam plot'
            beamAx = gca;
        subplot(4,3,[2,3,5,6])
            hold on
            grid on
            xlabel 'time [s]'
            ylabel 'displasement [mm]'
            title 'Laser Measurement'
            measurementAx = gca;
            xlim([0,simulationSettings.T]);
        subplot(4,3,[8,9,11,12])
            bodeAx = gca;
            
        % Plot bode
        bod = bodeplot(bodeAx,SS_modal(1,1));
        
        plotoptions = getoptions(bod);
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
        plotOptions.XLim = {[0,wmax/(2*pi)]};
      
        setoptions(bod,plotoptions);
            
        simPlots = [];
%         measurement = zeros(2,length(t));        
        
        if simulationSettings.animate == true
            Nsteps = length(t);
        else
            Nsteps = 1;
        end

        % Run animation loop! (Only once if animate is turned off)
        for i = 1:Nsteps
            tstart = tic;
            d = Phi*z(1:Nmodes,i);

            % update bodies
            for j = 1:length(bodies)
                bodies(j).update(d);
            end
            
            % Delete old plot
            if ~isempty(simPlots)
                delete(simPlots)
            end

            % Plot beam
            simPlots = plotModel(bodies,plotSettings,simulationSettings,beamAx); 
            timeText = text(beamAx, 0,beam.L*1.1,['Time: ',num2str(round(t(i),1)),'/',num2str(t(end)),' s'],...
                                                        'HorizontalAlignment','center');
            simPlots = [simPlots, timeText];
            
            % Plot laser
            if plotSettings.sensor == true
                laser = plot(beamAx,[beamAx.XLim(1) measurement(i)-beam.h/2],[1,1]*simulationSettings.measurementHeight*beam.L,'r','lineWidth',2);
                simPlots = [simPlots, laser];
            end
            
            % Plot laser measurement
            if plotSettings.sensorPlot == true
                if simulationSettings.animate == true
                    if i > 1
                        xs = [t(i-1),t(i)];
                        ys = [measurement(i-1),measurement(i)];
                        if plotSettings.base == true
                            yb = [baseLocation(i-1),baseLocation(i)];
                        end
                    else
                        xs = t(i);
                        ys = measurement(i);
                        yb = baseLocation(i);
                    end
                else
                    xs = t;
                    ys = measurement;
                    yb = baseLocation;
                end
                    if plotSettings.base == true
                        basePlot = plot(measurementAx,xs,yb,'r');
                    end
                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
            end
            
            % Plot input and output as dots
            forceAlpha = simulationSettings.forcePoint(1,4);
            inputNode = simulationSettings.forcePoint(:,[2,3])';
            inputNode = inputNode(:,1)+(inputNode(:,2)-inputNode(:,1))*forceAlpha;

            inputN = plot(beamAx,inputNode(1),inputNode(2),'g.','MarkerSize',12,'color',[0.4660 0.6740 0.1880]);
            simPlots = [simPlots, inputN];
            
            measureAlpha = simulationSettings.measurePoint(1,4);
            measureNode = simulationSettings.measurePoint(:,[2,3])';
            measureNode = measureNode(:,1)+(measureNode(:,2)-measureNode(:,1))*measureAlpha;
            
            outputN = plot(beamAx,measureNode(1),measureNode(2),'r.','MarkerSize',12);
            simPlots = [simPlots, outputN];
            
            drawnow;
            
            % Measure elapsed time and delay
            elapsed = toc(tstart);
            
            if simulationSettings.animate == true
                pause(max(0.01,0));
            end
            if simulationSettings.pauseStart == true && i == 1 && simulationSettings.animate == true
                disp('PAUSED, Press any key to continue');
                pause;
            end
        end 
    end
end

%% Mode shape analysis
if plotSettings.modalAnalysis == true
disp('Mode analysis plotting...')
% Plot all mode shapes
    if Nmodes <= 3
        i = 1;
        j = Nmodes;
    elseif Nmodes > 3 && Nmodes <= 9
        i = ceil(Nmodes/3);
        j = 3;
    elseif Nmodes > 9 && Nmodes <= 12
        i = 4;
        j = 3;
    elseif Nmodes > 12 && Nmodes <= 16
        i = 4;
        j = 4;
    elseif Nmodes > 16 && Nmodes <= 20
        i = 4;
        j = 5;
    else
        warning('Too many mode shapes to plot!')
    end

    figure()
    sgtitle(['First ',num2str(Nmodes),' eigenmodes']);
    hold on
    for k = 1:min(Nmodes,20)
        subplot(i,j,k)
            title(['Mode ',num2str(k),': ',num2str(round(modal.NaturalFrequencies(k)/(2*pi),2)),' Hz'])
            hold on
            xlabel('m')
            ylabel('m')
            grid on
            axis equal
            ax = gca;

            % Amplification factor
            A = plotSettings.A;
            % Update Beam
            updateBeam(elements,nodes,faces,A*(Phi(:,k)/norm(Phi(:,k))));
            % Plot Beam
            plotBeam(nodes,elements,faces,plotSettings,simulationSettings,ax);
    end
end
disp('Done!----------------------------------------------------------------')

%% Functions
% odefun for ode45 etc..
function qd = sysFun(t,q,A,B,simulationSettings)
    % Generic input!
    switch simulationSettings.Input 
        case {'force'}
            r = 0;
            if t > simulationSettings.stepTime
                F = 1;
                Q = 0;%1e10;
            else
                F = 0;
                Q = 0;
            end
        case {'disp'}
            F = 0;
            if t > simulationSettings.stepTime
                r = 1;
                Q = 0;
            else
                r = 0;
                Q = 0;
            end
    end        
    % Evaluate system
    qd = A*q+B*[F,r,Q*eye(1,simulationSettings.Nactuators)]';
end

% Function to update all nodes, links and faces.
function [links,nodes] =  updateBeam(links,nodes,faces,q)
    numNodes = length(nodes);
    % Nodes
    for i = 1:numNodes
        nodes{i}.update(q);
    end
    % Links
    for i = 1:length(links)
        links{i}.update(nodes);
    end
    % Faces
    for i = 1:length(faces)
        faces{i}.update(nodes);
    end
end   

% Function to find nearest nodes when height form the left is given. 
function interpNodes = interpolate(bodies,y)
    nodesxy = [bodies(1).nodeLocations];
    nodesxy = [1:length(nodesxy);
                nodesxy];
    
    % Nodes bove and below measurement height
    higher = nodesxy(:,nodesxy(3,:)>=y);
    lower =  nodesxy(:,nodesxy(3,:)<y);
        
    % sort on distance from measurement height
    higher = sortrows(higher',3,'ascend');
    lower = sortrows(lower',3,'descend');
    
    % Force to look on left edge of beam
    higher = higher(higher(:,2)==-bodies(1).h/2,:);
    lower = lower(lower(:,2)==-bodies(1).h/2,:);
    
    interpNodes = zeros(2,4);
    interpNodes(1,1:3) = lower(1,:);   % closest node below
    interpNodes(2,1:3) = higher(1,:);    % closest node above
    
    if interpNodes(1,3) == interpNodes(2,3)
        alpha = 0.5;
        
    else
        alpha = (y-interpNodes(1,3))/(interpNodes(2,3)-interpNodes(1,3));
    end
    
    interpNodes(:,4) =  [alpha;
                         alpha];
end
              
              
              