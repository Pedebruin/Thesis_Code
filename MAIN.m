% This is the main file.
clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

disp('Setting up...')
% Create template objects
beam = body('Beam');
sensor = body('Sensor');
actuator = body('Actuator');

%% Parameters & Settings
% beam parameters
beam.L = 370e-3;                    % m
beam.h = 2e-3;                      % m
beam.b = 40e-3;                     % m
beam.N = 65; %beam.L/beam.h;                   % Number of Vertical Nodes (minimum 2) (Accurate from around 65)
beam.E = 70E9;                      % E modulus
beam.G = 25.5E9;                    % Shear Modulus
beam.mu = 0.334;                    % Poisson
beam.rho = 2710;                    % Mass density
beam.alphaC = 0.001599;             % For proportional damping
beam.betaC = 0.001595;              % For proportional damping
beam.zeta = 0.01;                   % For modal damping

% sensor parameters PI P-876 DuraAct
simulationSettings.Nsensors = 2;    % number of sensor patches
sensor.L = 50e-3;                   % m
sensor.h = 500e-6; %1;                % m
sensor.b = 40e-3;                   % m
sensor.N = 2; %sensor.L/beam.h;
sensor.E = 5.2e10;
sensor.G = 0.775e9;                 % Shear modulus
sensor.mu = 0.334;
sensor.rho = 7800;
sensor.alphaC = 0.001599;
sensor.betaC = 0.001595;
sensor.zeta = 0.01;

sensor.d31 = -180e-12;
sensor.d33 = 400e-12;
sensor.s11 = 16.1e-12;
sensor.eta33 = 1750*sensor.eta0;

% actuator parameters PI P-876 DuraAct
simulationSettings.Nactuators = 4;  % number of actuators
actuator.L = 50e-3;                 % m
actuator.h = 500e-6; %1               % m
actuator.b = 40e-3;                 % m
actuator.N = 2; % actuator.L/beam.h;
actuator.E = 5.2e10;
actuator.G = 0.775e9;               % Shear modulus
actuator.mu = 0.334;
actuator.rho = 7800;
actuator.alphaC = 0.001599;
actuator.betaC = 0.001595;
actuator.zeta = 0.01;

actuator.d31 = -180e-12;
actuator.d33 = 400e-12;
actuator.s11 = 16.1e-12;
actuator.eta33 = 1750*actuator.eta0;

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
simulationSettings.forceHeight = 0.85;             % []*L Height of applied force
simulationSettings.dampingModel = 'modal';      % 'Modal','proportional' or 'none'
simulationSettings.Input = 'force';             % 'force' for force input, or 'disp' for disp input   
    simulationSettings.Kbase = 1e10;            % only in 'disp' mode

% Plot settings
% Single mode plot settings
plotSettings.modalAnalysis = false;          % Single plot with a specific mode
    plotSettings.modes = 9;
    plotSettings.A = 0.1;
    
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

    % For piezo augmentation, the states need to be augmented-> q = [z;zd;Q;Qd]
    patches = [sensors, actuators];

    % Piezo parameters
    [e33s,eta33s] = sensor.setupPiezo(bodies);
    [e33a,eta33a] = actuator.setupPiezo(bodies);
    ds = [sensor.d33, sensor.d31, 0;       % Piezo coupling matrix for PZT in 2D (NEED TO CHECK!)(Depends on direction of polarisation)
        0, 0, 0];
    es = ds/sensor.s11;
    etas = eye(2)*eta33s;

    da = [actuator.d33, actuator.d31, 0;       % Piezo coupling matrix for PZT in 2D (NEED TO CHECK!)(Depends on direction of polarisation)
        0, 0, 0];
    ea = da/actuator.s11;
    etaa = eye(2)*eta33a;

    T_s = sparse(zeros(simulationSettings.Nsensors,numNodes*2));
    Kpu_s = sparse(zeros(numNodes,numNodes*2));
    Kpp_s = sparse(zeros(numNodes,numNodes));
    
    T_a = sparse(zeros(simulationSettings.Nactuators,numNodes*2));
    Kpu_a = sparse(zeros(numNodes,numNodes*2));
    Kpp_a = sparse(zeros(numNodes,numNodes));
    
    for i = 1:Np                        % For every patch
        elements = [patches(i).elements];
        
        switch patches(i).name([1:3])   % is sensor or actuator patch?
            case {'Sen'}
                e = es;
                eta = etas;
                h = sensor.h;
            case {'Act'}
                e = ea;
                eta = etaa;
                h = actuator.h;
        end

        for j = 1:length(elements)      % For every element in the current patch
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

            Kpu = A*Bphi'*e*Bu;   % Need to look at this ds to derive for 2D!! ('omitted now)
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

            % Put matrix in larger matrix! (For sensors and actuators)
            switch patches(i).name([1:3])
                case {'Sen'}
                    Kpu_s = Kpu_s+Kpu_e;
                    Kpp_s = Kpp_s+Kpp_e;
                    T_s(i,[n1,n2,n3,n1+numNodes,n2+numNodes,n3+numNodes]) = ones(1,6);
                case {'Act'}
                    Kpu_a = Kpu_a+Kpu_e;
                    Kpp_a = Kpp_a+Kpp_e;
                    T_a(i-Nsensors,[n1,n2,n3,n1+numNodes,n2+numNodes,n3+numNodes]) = ones(1,6);
            end
        end
    end

    % Add all potentials on the nodes for the actuators and the sensors to
    % obtain the potential state of each patch! (Thomas, 2012)
    Kc = Kpu_s+Kpu_a;
    Ke = Kpp_s+Kpp_a;
    
    % actuator input and output vectors --->>> q = [z;zd;Q;Qd]
    Bact = T_a';    
    Csens = T_s;    

    numSensorNodes = sum(T_s,2);
    numActuatorNodes = sum(T_a,2);
    
    A = [Amech,                                 [zeros(Nmodes,numNodes),    zeros(Nmodes,numNodes);
                                                 -Phi'*Kc',                 zeros(Nmodes,numNodes)];   
        zeros(numNodes,Nmodes),   zeros(numNodes,Nmodes),   zeros(numNodes,numNodes),       eye(numNodes);
        Kc*Phi,            zeros(numNodes,Nmodes),   -Ke,                 zeros(numNodes)];
    B = [Bmech,               zeros(Nmodes*2,Nactuators);
        zeros(numNodes*2,size(Bmech,2)),    Bact];                           %Q
    C = [Cmech,                     zeros(1,numNodes*2);                   %[Laser;Sensors;Modal states];
         zeros(Nsensors,Nmodes*2),  Csens;
         eye(Nmodes),zeros(Nmodes),zeros(Nmodes,numNodes*2)];
    D = zeros(size(Cmech,1)+Nsensors+Nmodes,size(Bmech,2)+Nactuators);
    
    for i = 1:Nactuators
        B(:,i+1) = -B(:,i)/sum(B(:,i));% *sigma;    % TODO -> CHeck sigma!
    end
%     for i = 1:Nsensors
%         C(i,:) = 
%     end

    SS_tot = ss(full(A),full(B),full(C),D);
    disp('Constructing minimal realisation...');
    SS_tot = minreal(SS_tot); % Minimal realisation, but physical states are lost!
    warning('Physical state interpretation lost!');
else
    SS_tot = SS_modalr;
end

%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if simulationSettings.simulate ==  true
    disp('Simulating...')
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    
    % Time simulation!
    Asim = SS_tot.A;
    Bsim = SS_tot.B;
    q0 = zeros(size(Asim,1),1);
    [t,y] = ode45(@(t,z) sysFun(t,z,Asim,Bsim,simulationSettings),tvec,q0);
    z = y';
    y = SS_tot.C*z;
    
    d = Phi*y(Nsensors+2:end,:);
    baseLocation = d(baseNodes(2),:);
    measurement = y(1,:);


   switch simulationSettings.Input
        case {'force'}
            plotSettings.base = false;
        case {'disp'}
    end
    
    if simulationSettings.plot == true
        disp('Plotting simulation...')
        figure()
        %sgtitle(['Beam with ',simulationSettings.Input,' input'])
        subplot(7,3,[1,4,7,10])
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
        subplot(7,3,[2,3,5,6])
            hold on
            grid on
            ylabel 'displasement [mm]'
            title 'Laser Measurement'
            measurementAx = gca;
            xlim([0,simulationSettings.T]);
        subplot(7,3,[8,9,11,12]);
            hold on
            grid on
            xlabel 'time [s]'
            ylabel 'Modal contribution'
            title 'Modal state evolution'   
            stateAx = gca;
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
        plotOptions.XLim = {[0,wmax/(2*pi)]};
        
        switch simulationSettings.Input
            case {'force'}
                bod = bodeplot(bodeAx,minreal(SS_modalr(1,1)),plotoptions);
            case {'disp'}
                bod = bodeplot(bodeAx,minreal(SS_modalr(1,2)),plotoptions);
        end
        
        
        % labels = {'1st','2nd','3rd','4th','5th','6th','7th'};
        % legend(stateAx,labels(1:min(Nmodes,3)));
        
        % Plot analytical omega1
        xline(bodeAx,analyticalOmega1)
      
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
                d = d(:,i);
            else
                d = d(:,end);
            end

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
                        yst = [z(1:Nmodes,i-1),z(1:Nmodes,i)];
                        if plotSettings.base == true
                            yb = [baseLocation(i-1),baseLocation(i)];
                        end
                    else
                        xs = t(i);
                        ys = measurement(i);
                        yb = baseLocation(i);
                        yst = y(Nsensors+2:end,i);
                    end
                else
                    xs = t;
                    ys = measurement;
                    yb = baseLocation;
                    yst = y(Nsensors+2:end,:);
                end
                    if plotSettings.base == true
                        basePlot = plot(measurementAx,xs,yb,'r');
                    end
                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
                    statePlot = plot(stateAx,xs,yst);
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
    
if strcmp(simulationSettings.Input,'disp')
    warning('Modes plotted are before the addition of the stiff spring!, have to analyse this later')
end

disp('Mode analysis plotting...')
Nmodesp = plotSettings.modes;

% Plot all mode shapes
    if Nmodesp <= 3
        i = 1;
        j = Nmodesp;
    elseif Nmodesp > 3 && Nmodesp <= 9
        i = ceil(Nmodesp/3);
        j = 3;
    elseif Nmodesp > 9 && Nmodesp <= 12
        i = 4;
        j = 3;
    elseif Nmodesp > 12 && Nmodesp <= 16
        i = 4;
        j = 4;
    elseif Nmodesp > 16 && Nmodesp <= 20
        i = 4;
        j = 5;
    else
        warning('Too many mode shapes to plot!')
    end

    figure()
    sgtitle(['First ',num2str(Nmodesp),' eigenmodes']);
    hold on
    for k = 1:min(Nmodesp,20)
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

            % update bodies
            for l = 1:length(bodies)
                bodies(l).update(A*Phi(:,k));
            end
            
            % Plot Beam
            plotModel(bodies,plotSettings,simulationSettings,ax);
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
                F = 0;
                Q = 1e2;
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
              
              
              