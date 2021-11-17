%{
This is a 1D Euler bernoulli simulation with piëzo patches

%}
clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

disp('Setting up...')

%% Parameters & Settings
L = 370e-3;     % Total length
b = 5e-3;      % Total width

% Model settings
modelSettings.numEl = 5;
modelSettings.sElements = [];   % Write a smarter function for this later!
modelSettings.dispInput = false;
modelSettings.Nmodes = 4;
modelSettings.zeta = 0.01;
modelSettings.measurementHeight = 0.5;        

% plotSettings
plotSettings.plotNodes = true;
    plotSettings.nodeNumbers = true;
plotSettings.elementNumbers = true;
plotSettings.sensor = true;
plotSettings.sensorPlot = true;

% simulationSettings
simulationSettings.simulate = true;
simulationSettings.stepTime = 0.1;
simulationSettings.dt = 0.001;
simulationSettings.T = 1;

simulationSettings.animate = false;
simulationSettings.plot = true;
simulationSettings.pauseStart = false;


% beam element parameters
Beam = element('Beam');
Beam.L = L/modelSettings.numEl;
Beam.h = 1e-3;                      % m
Beam.b = b;                         % m
Beam.E = 70E9;                      % E modulus
Beam.mu = 0.334;                    % Poisson
Beam.rho = 2710;                    % Mass density

Beam.N = 10;                        % Number of interpolation points for plotting

% Smart beam parameters :PI P-876 DuraAct
sBeam = copy(Beam); sBeam.name = 'sBeam';
sBeam.ph = 0.5e-3;                  % m
sBeam.pb = b;                       % m
sBeam.pE = 5.2e10;
sBeam.pmu = 0.334;
sBeam.prho = 7800;

%% Matrix setup 
% Assembly
numEl = modelSettings.numEl;
numNodes = numEl+1;
Nmodes = modelSettings.Nmodes;

K = zeros(numNodes*2);
M = zeros(numNodes*2);
elements = element.empty(numEl,0);

% Loops over all elements, constructs their elemental K and M matrices and
% stores them in element objects aswell. Also constructs the system K and M
% matrices. 
for i = 1:numEl
    n1 = i;     % Starting node
    n2 = i+1;   % Ending node
  
    if sum(ismember(modelSettings.sElements,i)) > 0 % Is smart element?
        El = copy(sBeam);
        El.sBeam = true;
    else                                            % Or normal beam element?
        El = copy(Beam);
    end
    
        Ib = El.b*El.h^3/12;
        I = El.b*El.ph^3/12 + El.b*El.ph*(El.ph+El.h)*2/4;
        
        EI = El.E*Ib+2*El.pE*I;                     % Equivalent EI and Rho*A (for normal element, extra term is 0)
        rhoA = El.b*(El.rho*El.h+2*El.prho*El.ph);

        Ke = EI/El.L^3*    [12, 6*El.L, -12, 6*El.L;
                            6*El.L, 4*El.L^2, -6*El.L, 2*El.L^2;
                            -12, -6*El.L, 12, -6*El.L;
                            6*El.L, 2*El.L^2, -6*El.L, 4*El.L^2];

        Me = rhoA*El.L/420* [156, 22*El.L, 54, -13*El.L;
                            22*El.L, 4*El.L^2, 13*El.L, -3*El.L^2;
                            54, 13*El.L, 156, -22*El.L;
                            -13*El.L, -3*El.L^2, -22*El.L, 4*El.L^2];    
    
    % Apply boundary condition to first node
    if n1 == 1
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
        end
    end
    
    % Assembly. Since everything is nicely connected in sequence, the
    % elemenental matrices stay nicely together.
    K(n1*2-1:n2*2,n1*2-1:n2*2) = K(n1*2-1:n2*2,n1*2-1:n2*2) + Ke;   % K
    M(n1*2-1:n2*2,n1*2-1:n2*2) = M(n1*2-1:n2*2,n1*2-1:n2*2) + Me;   % M
    
    %% Put information in elements
    x1 = 0;
    x2 = 0;
    y1 = (i-1)*Beam.L;
    y2 = i*Beam.L;
    t1 = 0;
    t2 = 0;
    initPos = [x1, x2;
                y1, y2;
                t1, t2];
    
    elements(i) = El.assign(i,[n1,n2],initPos,Ke,Me);
end

%% 
if modelSettings.Nmodes <= modelSettings.numEl
    [Phi,omega2] = eigs(K,M,modelSettings.Nmodes+2,'smallestabs');
else
    error('Number of elements is smaller then the number of modes requested. Modal decomposition will not work!')
end

% Snip off weird rigid body modes (Probably due to the setting to identity
% of fixed DOFs, and not deleting them). 
Phi = Phi(:,3:end);
omega2 = omega2(3:end,3:end);

% Check eigenfrequencies
f = diag(sqrt(omega2)/(2*pi));

I = (Beam.h^3)/12;
S = Beam.h;         % Need to check this formulation still!
lambda = 1.87510407;
analyticalf1 = lambda^2/(2*pi*L^2)*sqrt(Beam.E*I/(Beam.rho*S));

if abs(f(1) - analyticalf1) > 0.001
    warning('Eigenfrequencies do not correspond')
end

%% Dynamical matrices setup

Phi = Phi/(Phi'*M*Phi);     % Normalise w.r.t. mass matrix

% Sensor equation according to K. Aktas
H = 1e6;
z = sBeam.h/2+sBeam.ph;
d31 = -180e-12;
s11 = 16.1e-12;
e31 = d31/s11;
w = sBeam.pb;
S = H*z*e31*w*[0,-1,0,1];

% Actuator equation according to K. Aktas
Ep = sBeam.pE;
d31 =  -180e-12;
w = sBeam.pb;
zbar = (sBeam.ph+sBeam.h)/2;
G = Ep*d31*w*zbar*[-1,0,1,0]';

% Damping matrix
Cmodal = 2*modelSettings.zeta*sqrt(omega2);

% External input B matrix (in d coordinates)
Bext = zeros(numNodes*2,1);
Bext(end-1) = 1;                          % Force at last node in x direction 

% Measurement matrix C (in d coordinates)
height = modelSettings.measurementHeight*L;
pos = unique([elements.pos]');
diff = pos-height;
lowerNodes = find(diff<0);
interpNodes = [lowerNodes(end),lowerNodes(end)+1];

for i = 1:numEl
    if interpNodes == elements(i).neighbours
        interpEl = i;
    end
end

alpha = height - pos(interpNodes(1));
N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
Cmeasurement = zeros(1,numNodes*2);
Cmeasurement(1,interpNodes(1)*2-1:interpNodes(1)*2+2) = N;

% Voltage input B matrix (in d coordinates)
nsElements = length(modelSettings.sElements);
Bg = zeros(numNodes*2,nsElements);
Cs = zeros(nsElements,numNodes*2);
for i = 1:nsElements
    el = modelSettings.sElements(i);
    n1 = elements(el).neighbours(1);
    
    Bg(n1*2-1:n1*2+2,i) = G;
    Cs(i,n1*2-1:n1*2+2) = S;
end

A = [zeros(Nmodes),eye(Nmodes);
    -omega2,-Cmodal];
B = [zeros(Nmodes,nsElements+1);
    Phi'*Bext,Phi'*Bg];
C = [Cmeasurement*Phi, zeros(1,Nmodes);
    Cs*Phi, zeros(nsElements,Nmodes);
    eye(Nmodes),zeros(Nmodes)];                           % Also output full state (for plotting)

sys = ss(A,B,C,[]);


%% Time simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if simulationSettings.simulate ==  true
    disp('Simulating...')
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    q0sim = zeros(size(A,1),1);
    
        [t,y] = ode45(@(t,z) sysFun(t,z,A,B,modelSettings,simulationSettings),tvec,q0sim);
        z = y';
        y = C*z;
    
    dmat = Phi*z(1:Nmodes,:);
    measurement = y(1,:);
    
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
  
        bod = bodeplot(bodeAx,sys(1,1),plotoptions);

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
            
            % Plot laser measurement
            if plotSettings.sensorPlot == true
                if simulationSettings.animate == true
                    if i > 1
                        xs = [t(i-1),t(i)];
                        ys = [measurement(i-1),measurement(i)];
                        yst = [z(1:Nmodes,i-1),z(1:Nmodes,i)];
                    else
                        xs = t(i);
                        ys = measurement(i);
                        yst = z(1:Nmodes,i);
                    end
                else
                    xs = t;
                    ys = measurement;
                    yst = z(1:Nmodes:end,:);
                end
                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
                    statePlot = plot(stateAx,xs,yst);
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

d = Phi(:,1);
d = d/max(abs(d));

figure()
hold on
axis equal
for i = 1:numEl
    elements(i).update(d);
    elements(i).show([],plotSettings);
end

disp('Done!')

% odefun for ode45 etc..
function qd = sysFun(t,q,A,B,modelSettings,simulationSettings)
    % Generic input!
    nsElements = length(modelSettings.sElements);
    if t > simulationSettings.stepTime
        F = 0.1;
        V = 0;%100;    % Input voltage on the patch (I think)
    else
        F = 0;
        V = 0;
    end

    % Evaluate system
    U = [F,V*eye(1,nsElements)]';
    qd = A*q+B*U;

    if sum(isnan(q))>0 || sum(isnan(qd)) > 0
        error("NaN's found in state vector")
    end        
end 