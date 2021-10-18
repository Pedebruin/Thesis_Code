clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

%% Parameters & Settings
% Parameters
L = 370;                    % mm
h = 2;                      % mm
b = 40;                     % mm

N = L/h;               % Number of Vertical Nodes (minimum 2)

E = 70E9;                   % E modulus
mu = 0.3;                   % Poisson
rho = 2700;                 % Mass density
alpha = 0.001599;           % For proportional damping
beta = 0.001595;            % For proportional damping
zeta = 0.01;                % For modal damping

wmax = 10000;               % Maximum frequency to include in modal decomposition

% Simulation settings
simulationSettings.simulate = true;
    simulationSettings.plot = true;
    simulationSettings.animate = false;
        simulationSettings.pauseStart = false;
    simulationSettings.stepTime = 1;            % s
    simulationSettings.dampingModel = 'modal';  % 'Modal','proportional' or 'none'
    simulationSettings.Kbase = 1e10;
    simulationSettings.T = 15;          % s
    simulationSettings.dt = 0.01;       % s 
    
    simulationSettings.measurementHeight = 0.75;    % x*L
        
% Single mode plot settings
plotSettings.modalAnalysis = false;          % Single plot with a specific mode
    plotSettings.modeNr = 2;
    plotSettings.A = 1e4;

% Plot settings
plotSettings.realMesh = false;                  % Make an extra plot with the original Mesh
plotSettings.color = 'none';                    % 'disp' for displacement or 'deform' for deformation or 'none'
    plotSettings.fill = true;                   % Fill the beam with patch
plotSettings.nodes = false;                     % Plot the nodes 
    plotSettings.nodeColor = 'k';               % Node color
    plotSettings.nodeNumbers = false;           % Node numbers
plotSettings.links = true;                     % Plot the links (edges)
    plotSettings.linkNumbers = false;           % link numbers
plotSettings.sensor = true;                     % plot the sensor
plotSettings.sensorPlot = true;                 % plot the sensor output


%% build FEM model of the beam
[FEM,nodes,links,faces,beam] = buildBeam(L,h,N,E,mu,rho,plotSettings);
numNodes = length(nodes);
numLinks = length(links);
numFaces = length(faces);

K = FEM.K;
M = FEM.M;

% Get modal results
wmin = -0.1;

% Modal decomposition! P395
modal = solve(beam,'FrequencyRange',[wmin,wmax]);   % Solve
Nmodes = length(modal.NaturalFrequencies);

I = h^3/12;
analyticalOmega1 = 3.516*sqrt(E*I/(L^4*(rho*h)))/(2*pi);

% Find node to sense!
nodesxy = beam.Mesh.Nodes;
measurey = simulationSettings.measurementHeight*L;
temp = [nodesxy(1,:);
        abs(nodesxy(2,:)-measurey);
        1:numNodes;];
temp = sortrows(temp',2);
temp = temp(temp(:,1)<0,:);
interpNodes = temp([1,2],3);
interpNodes = sortrows([interpNodes,nodesxy(2,interpNodes)'],2);
alpha = (measurey-interpNodes(1,2))/(interpNodes(2,2)-interpNodes(1,2));

simulationSettings.interpPoint = [interpNodes,[alpha;
                                        alpha]];
                                                                     
% Setup Model
% Input
Bd = zeros(numNodes*2,1);                                   % External Input matrix in d space
Bd(1) = 1;                                                  % Force on node 1 in x;                                                 
Bd(2) = 1;                                                  % Force on node 2 in x;

% Output
Cd = zeros(1,numNodes*2);                                   % Measurement matrix in d space
Cd(interpNodes(:,1)) = [1-alpha;alpha];                     % Measure at selected point

%% modal system modelling
Phi = [modal.ModeShapes.ux;                         % Phi
    modal.ModeShapes.uy];
Phi = Phi/(Phi'*M*Phi);

omegaSq = Phi'*K*Phi;                                           % Omega^2

% Damping matrix C_Phi
switch simulationSettings.dampingModel                          % Damping matic C_phi
    case {'modal'}
        C_phi = diag(2*zeta*modal.NaturalFrequencies);
    case {'proportional'}
        C_phi = alpha*eye(Nmodes)+beta*omegaSq;
    case {'none'}
        C_phi = zeros(Nmodes);
end
 
Rphi = Phi'*Bd;                                             % Modal input matrix
Cphi = Cd*Phi;                                              % Modal measurement

% Assemble new system for modal
A_modal = [zeros(Nmodes), eye(Nmodes);                            % Dynamics
     -omegaSq,-C_phi];
B_modal = [zeros(Nmodes,1);Rphi];                                 % External input
C_modal = [Cphi,zeros(1,Nmodes)];                                 % Measurement
D_modal = 0;                                                % Assume no direct feedthrough

SS_force = ss(A_modal,B_modal,C_modal,D_modal);

% Apply stiff spring to base   
Ck = zeros(1,numNodes*2);               % Ficticious measurement matrix
Ck(1) = 1;                              % measure displacement of node 1
K_spring = simulationSettings.Kbase*[Ck*Phi,zeros(1,Nmodes)];    % get stiffness matrix

A_modal = A_modal-B_modal*K_spring;   % Apply feedback to simulate spring

% Get system! The system now has a state vector q = [z;zdot]; and an input
% r which is the reference location of node 1. 
SS_modal = ss(A_modal,B_modal,C_modal,D_modal);
gain = dcgain(SS_modal);
B_modal = B_modal/gain;
SS_modal.B = B_modal;

%% Time simulation
if simulationSettings.simulate ==  true
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    
    % Time simulation!
    q0 = zeros(Nmodes*2,1);
    [t,y] = ode45(@(t,z) sysFun(t,z,A_modal,B_modal,simulationSettings),tvec,q0);
    z = y';
    output = C_modal*z;
    baseLocation = [Ck*Phi,zeros(1,Nmodes)]*z;
    
    if simulationSettings.plot == true
        figure()
        subplot(4,3,[1,4,7,10])
            hold on
            axis equal
            grid on
            xlabel 'mm'
            ylabel 'mm'
            xlim([-L/6,L/6])
            ylim([0,1.2*L])
            title 'Beam plot'
            beamAx = gca;
        subplot(4,3,[2,3,5,6])
            hold on
            grid on
            xlabel 'time [s]'
            title 'Laser Measurement'
            measurementAx = gca;
            xlim([0,simulationSettings.T]);
        subplot(4,3,[8,9,11,12])
            bodeAx = gca;
            
        % Plot bode
        h = bodeplot(bodeAx,SS_modal);
        
        plotoptions = getoptions(h);
        plotoptions.Title.String = '';
        plotoptions.Title.Interpreter = 'latex';
        plotoptions.XLabel.Interpreter = 'latex';
        plotoptions.YLabel.Interpreter = 'latex';
        plotoptions.XLabel.FontSize = 9;
        plotoptions.YLabel.FontSize = 9;
        plotoptions.FreqUnits = 'Hz';
      
        setoptions(h,plotoptions);
            
        simPlots = [];
        measurement = zeros(2,length(t));        
        
        if simulationSettings.animate == true
            Nsteps = length(t);
        else
            Nsteps = 1;
        end
        for i = 1:Nsteps
            tstart = tic;
            d = Phi*z(1:Nmodes,i);

            % Update Beam
            updateBeam(links,nodes,faces,d);
            
            % Delete old plot
            if ~isempty(simPlots)
                delete(simPlots)
            end

            % Plot beam
            simPlots = plotBeam(nodes,links,faces,plotSettings,simulationSettings,beamAx); 
            timeText = text(beamAx, 0,L*1.1,['Time: ',num2str(round(t(i),1)),'/',num2str(t(end)),' s'],...
                                                        'HorizontalAlignment','center');
            simPlots = [simPlots, timeText];
            
            % Plot laser measurement
            if plotSettings.sensor == true
                laser = plot(beamAx,[beamAx.XLim(1) output(i)],[1,1]*simulationSettings.measurementHeight*L,'r','lineWidth',2);
                simPlots = [simPlots, laser];
            end
            
            if plotSettings.sensorPlot == true
                if simulationSettings.animate == true
                    if i > 1
                        xs = [t(i-1),t(i)];
                        ys = [output(i-1),output(i)];
                    else
                        xs = t(i);
                        ys = output(i);
                    end
                else
                    xs = t;
                    ys = output;
                end

                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
            end
            
            drawnow;
            
            % Measure elapsed time and delay
            elapsed = toc(tstart);
            
            if simulationSettings.animate == true
                pause(max(simulationSettings.dt-elapsed,0));
            end
            if simulationSettings.pauseStart == true && i == 1 && simulationSettings.animate == true
                disp('PAUSED, Press any key to continue');
                pause;
            end
        end 
    end
end


%% Single mode plot
if plotSettings.modalAnalysis == true
    nr = plotSettings.modeNr;
    A = plotSettings.A;

    % Update Beam
    updateBeam(links,nodes,faces,A*Phi(:,nr));

    % Plot beam
    figure()
    hold on
    axis equal
    xlabel 'mm'
    ylabel 'mm'
    grid on
    title 'Beam plot'
    ax = gca;
    
    plotBeam(nodes,links,faces,plotSettings,simulationSettings,ax);
end


%% Functions
% odefun for ode45 etc..
function qd = sysFun(t,q,A,B,simulationSettings)
    % Apply Spring to base. with reference r
    if t > simulationSettings.stepTime
        r = 1;
    else
        r = 0;
    end
    
    % Evaluate system
    qd = A*q+B*r;
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
              
              
              