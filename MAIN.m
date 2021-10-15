clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

%% Parameters & Settings
% Parameters
L = 300;    % mm
h = 2;     % mm

N = 10;%L/h;      % Number of Vertical Nodes (minimum 2)

E = 210E9;
mu = 0.3;
rho = 7800;
alpha = 10;
beta = 2;
zeta = 0.1;

simulationSettings.simulate = true;
    simulationSettings.damping = true;
    simulationSettings.dampingModel = 'modal';
    simulationSettings.T = 25;       % s
    simulationSettings.dt = 0.1;   % s 
    simulationSettings.plot = true;

% Settings
plotSettings.realMesh = false;
plotSettings.modalAnalysis = false;
    plotSettings.modeNr = 3;
    plotSettings.A = 1e4;
plotSettings.nodes = true;
plotSettings.nodeNumbers = false;
plotSettings.links = true;
plotSettings.linkNumbers = false;


% build FEM model of the beam
[FEM,nodes,links,beam] = buildBeam(L,h,N,E,mu,rho,plotSettings);
numNodes = length(nodes);
numLinks = length(links);

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
if strcmp(simulationSettings.dampingModel,'modal')              % Cphi
    Cphi = diag(2*zeta*modal.NaturalFrequencies);
elseif strcmp(simulationSettings.dampingModel,'proportional')
    Cphi = alpha*eye(Nmodes)+beta*omegaSq;              
end
Rext = zeros(numNodes*2,1);                                 % External Input
Rphi = Phi'*Rext;                                         


%% Time simulation
if simulationSettings.simulate ==  true
    z0 = zeros(Nmodes,1);
    z0(2) = 4e4;
    z0(1) = 2e4;
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
        xlim([0,1.2*L])
        ylim([-L/6,L/6])
        title 'Beam plot'
        ax = gca;
        simPlots = [];
        for i = 1:length(t)
            d = Phi*z(:,i);

            % Update Beam
            updateBeam(links,nodes,d);

            if ~isempty(simPlots)
                delete(simPlots)
            end

            % Plot beam
            simPlots = plotBeam(nodes,links,plotSettings,ax); 
            P = text(ax, L/3*2,L/6*2,['Time; ',num2str(t(i)),'/',num2str(t(end)),' s']);
            simPlots = [simPlots, P];
            pause(0.01);
        end 
    end
end

%% Single mode plot
if plotSettings.modalAnalysis == true
    nr = plotSettings.modeNr;
    A = plotSettings.A;

    % Update Beam
    updateBeam(links,nodes,A*Phi(:,nr));

    % Plot beam
    figure()
    hold on
    axis equal
    xlabel 'mm'
    ylabel 'mm'
    title 'Beam plot'
    ax = gca;
    
    plotBeam(nodes,links,plotSettings,ax);
end

%% Functions
% odefun for ode45 etc..
function yd = sysFun(t,y,A,B)
% Kept the odefun as basic as possible, for future simulations!
yd = A*y+B;
end

function [links,nodes] =  updateBeam(links,nodes,q)
    numNodes = length(nodes);
    % Update beam
    % Nodes
    for i = 1:numNodes
        nodes{i}.update(q);
    end
    % Links
    for i = 1:length(links)
        links{i}.update(nodes);
    end
end             
              
              
              