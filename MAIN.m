clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

%% Parameters & Settings
% Parameters
L = 300;    % mm
h = 2;     % mm

N = 10; %L/h;      % Number of Vertical Nodes (minimum 2)

E = 210E9;
mu = 0.3;
rho = 7800;
alpha = 0.001599;         % For modal damping
beta = 0.001595;           % For modal damping
zeta = 0.005;        % For proportional damping

% Simulation settings
simulationSettings.simulate = true;
    simulationSettings.damping = true;
    simulationSettings.dampingModel = 'proportional';
    simulationSettings.T = 15;       % s
    simulationSettings.dt = 0.1;   % s 
    simulationSettings.plot = true;

% Single mode plot settings
plotSettings.modalAnalysis = true;          % Single plot with a specific mode
    plotSettings.modeNr = 1;
    plotSettings.A = 1e4;

% Plot settings
plotSettings.realMesh = false;               % Make an extra plot with the original Mesh
plotSettings.color = 'disp';        % 'disp' for displacement or 'deform' for deformation or 'none'
    plotSettings.fill = true;                   % Fill the beam with patch
plotSettings.nodes = false;
    plotSettings.nodeColor = 'k';
    plotSettings.nodeNumbers = false;
plotSettings.links = false;
    plotSettings.linkNumbers = false;


% build FEM model of the beam
[FEM,nodes,links,faces,beam] = buildBeam(L,h,N,E,mu,rho,plotSettings);
numNodes = length(nodes);
numLinks = length(links);
numFaces = length(faces);

K = FEM.K;
M = FEM.M;
R = FEM.R;

% Get modal results
wmin = -0.1;
wmax = 100;

% Modal decomposition! P395
modal = solve(beam,'FrequencyRange',[wmin,wmax]);   % Solve
Nmodes = length(modal.NaturalFrequencies);
Phi = [modal.ModeShapes.ux;                         % Phi
    modal.ModeShapes.uy];
Phi = Phi/(Phi'*M*Phi);

% Obtain modal system matrices
omegaSq = Phi'*K*Phi;                                           % Omega^2
switch simulationSettings.dampingModel
    case {'modal'}
        Cphi = diag(2*zeta*modal.NaturalFrequencies);
    case {'proportional'}
        Cphi = alpha*eye(Nmodes)+beta*omegaSq;
end

Rext = zeros(numNodes*2,1);                                 % External Input
Rphi = Phi'*Rext;                                         


%% Time simulation
if simulationSettings.simulate ==  true
    z0 = zeros(Nmodes,1);
    z0(2) = 4e3;
    z0(1) = 2e3;
    tvec = 0:simulationSettings.dt:simulationSettings.T;

    % y -> [z;zdot];!!!!!
    A = [zeros(Nmodes), eye(Nmodes);
         -omegaSq,-Cphi];
    B = [Rphi;zeros(Nmodes,1)];
    y0 = [z0;zeros(Nmodes,1)];

    [t,y] = ode45(@(t,z) sysFun(t,z,A,B),tvec,y0);
    z = y(:,1:Nmodes)';
    
    if simulationSettings.plot == true
        figure()
        hold on
        axis equal
        grid on
        xlabel 'mm'
        ylabel 'mm'
        xlim([-L/6,L/6])
        ylim([0,1.2*L])
        title 'Beam plot'
        ax = gca;
        simPlots = [];
        for i = 1:length(t)
            tstart = tic;
            d = Phi*z(:,i);

            % Update Beam
            updateBeam(links,nodes,faces,d);

            if ~isempty(simPlots)
                delete(simPlots)
            end

            % Plot beam
            simPlots = plotBeam(nodes,links,faces,plotSettings,ax); 
            P = text(ax, 0,L*1.1,['Time: ',num2str(round(t(i),1)),'/',num2str(t(end)),' s'],...
                                                        'HorizontalAlignment','center');
            simPlots = [simPlots, P];
            
            elapsed = toc(tstart);
            pause(max(simulationSettings.dt-elapsed,0));
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
    
    plotBeam(nodes,links,faces,plotSettings,ax);
end

%% Functions
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
              
              
              