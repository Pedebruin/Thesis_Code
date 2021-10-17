clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

%% Parameters & Settings
% Parameters
L = 300;    % mm
h = 2;     % mm

N = 5; %L/h;      % Number of Vertical Nodes (minimum 2)

E = 210E9;
mu = 0.3;
rho = 7800;
alpha = 0.001599;           % For proportional damping
beta = 0.001595;            % For proportional damping
zeta = 0.01;                % For modal damping

wmax = 100;                 % Maximum frequency to include in modal decomposition

% Simulation settings
simulationSettings.simulate = true;
    simulationSettings.plot = true;
    simulationSettings.pauseStart = true;
    simulationSettings.damping = true;
    simulationSettings.dampingModel = 'modal';
    simulationSettings.T = 15;          % s
    simulationSettings.dt = 0.1;       % s 
    
    
    simulationSettings.measurementHeight = 0.75;
        
% Single mode plot settings
plotSettings.modalAnalysis = false;          % Single plot with a specific mode
    plotSettings.modeNr = 1;
    plotSettings.A = 1e4;

% Plot settings
plotSettings.realMesh = false;                  % Make an extra plot with the original Mesh
plotSettings.color = 'none';                    % 'disp' for displacement or 'deform' for deformation or 'none'
    plotSettings.fill = true;                   % Fill the beam with patch
plotSettings.nodes = false;                     % Plot the nodes 
    plotSettings.nodeColor = 'k';               % Node color
    plotSettings.nodeNumbers = false;           % Node numbers
plotSettings.links = false;                     % Plot the links (edges)
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
R = FEM.R;

% Get modal results
wmin = -0.1;

% Modal decomposition! P395
modal = solve(beam,'FrequencyRange',[wmin,wmax]);   % Solve
Nmodes = length(modal.NaturalFrequencies);


%% Find node to sense!
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
                                    
                                    
%% Setup Model

% Apply input and measurement
Bd = zeros(numNodes*2,1);                                   % External Input matrix in d space
Cd = zeros(1,numNodes*2);                                   % Measurement matrix in d space

Bd(3) = 1;                                                  % For now only force on node 3 
Cd(interpNodes(:,1)) = [1-alpha;alpha];                     % Measure at selected point
Phi = [modal.ModeShapes.ux;                         % Phi
    modal.ModeShapes.uy];
Phi = Phi/(Phi'*M*Phi);

% Obtain modal system matrices
omegaSq = Phi'*K*Phi;                                           % Omega^2
switch simulationSettings.dampingModel                          % Damping matic C_phi
    case {'modal'}
        C_phi = diag(2*zeta*modal.NaturalFrequencies);
    case {'proportional'}
        C_phi = alpha*eye(Nmodes)+beta*omegaSq;
end
 
Rphi = Phi'*Bd;                                             % Modal input matrix
Cphi = Cd*Phi;                                              % Modal measurement

% Assemble new system
A = [zeros(Nmodes), eye(Nmodes);                            % Dynamics
     -omegaSq,-C_phi];
B = [Rphi;zeros(Nmodes,1)];                                 % External input
C = [Cphi,zeros(1,Nmodes)];                                 % Measurement
D = 0;                                                     % Assume no direct feedthrough

% Get system
SSPhi = ss(A,B,C,D);

%% Time simulation
if simulationSettings.simulate ==  true
    z0 = zeros(Nmodes,1);
    z0(4) = 4e3;
    z0(2) = 4e3;
    z0(1) = 4e3;
    %z0 = 10*ones(Nmodes,1);
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    q0 = [z0;zeros(Nmodes,1)];  % q = [z;zdot]


    [t,y] = ode45(@(t,z) sysFun(t,z,A,B),tvec,q0);
    z = y';
    
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
            xlabel 'time [s]'
            title 'Laser Measurement'
            measurementAx = gca;
            xlim([0,simulationSettings.T]);
        subplot(4,3,[8,9,11,12])
            bodeAx = gca;
            
        % Plot bode
        h = bodeplot(bodeAx,SSPhi);
        
        plotoptions = getoptions(h);
        plotoptions.Title.String = '';
        plotoptions.Title.Interpreter = 'latex';
        plotoptions.XLabel.Interpreter = 'latex';
        plotoptions.YLabel.Interpreter = 'latex';
        plotoptions.XLabel.FontSize = 9;
        plotoptions.YLabel.FontSize = 9;
      
        setoptions(h,plotoptions);
            
        % Time simulation!
        simPlots = [];
        measurement = zeros(2,length(t));
        output = zeros(1,length(t));
        
        
        
        for i = 1:length(t)
            tstart = tic;
            d = Phi*z(1:Nmodes,i);

            % Update Beam
            updateBeam(links,nodes,faces,d);
            
            % Measure
            output(i) = C*z(:,i);
            
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
                if i > 1
                    xs = [t(i-1),t(i)];
                    ys = [output(i-1),output(i)];
                
                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
                end
            end
            
            drawnow;
            
            % Measure elapsed time and delay
            elapsed = toc(tstart);
            pause(max(simulationSettings.dt-elapsed,0));
            
            if simulationSettings.pauseStart == true && i == 1
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
function [x,y] = measurePoint(simulationSettings,nodes)
        interpNodes = simulationSettings.interpPoint;
        alpha = interpNodes(1,3);
        n1 = nodes{interpNodes(1,1)}.pos;
        n2 = nodes{interpNodes(2,1)}.pos;
        
        measurePoint = n1+alpha*(n2-n1);
        x = measurePoint(1);
        y = measurePoint(2);
end

% odefun for ode45 etc..
function yd = sysFun(t,y,A,B)
% Kept the odefun as basic as possible, for future simulations!
yd = A*y+B;
end

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
              
              
              