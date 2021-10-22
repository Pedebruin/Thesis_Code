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
beam.N = 3;%65; %L/h;                  % Number of Vertical Nodes (minimum 2) (Accurate from around 65)
beam.E = 70E9;                      % E modulus
beam.mu = 0.334;                    % Poisson
beam.rho = 2710;                    % Mass density
beam.alphaC = 0.001599;             % For proportional damping
beam.betaC = 0.001595;              % For proportional damping
beam.zeta = 0.01;                   % For modal damping

% sensor parameters
sensor.L = 30e-3;                   % m
sensor.h = 1e-3;                    % m
sensor.N = 4;
sensor.E = 70e9;
sensor.mu = 0.334;
sensor.rho = 2710;
sensor.alphaC = 0.001599;
sensor.betaC = 0.001595;
sensor.zeta = 0.01;

% actuator parameters
actuator.L = 30e-3;                   % m
actuator.h = 1e-3;                    % m
actuator.N = 2;
actuator.E = 70e9;
actuator.mu = 0.334;
actuator.rho = 2710;
actuator.alphaC = 0.001599;
actuator.betaC = 0.001595;
actuator.zeta = 0.01;

wmax = 1e4*(2*pi);                    % Maximum frequency to include in modal decomposition

% Simulation settings
simulationSettings.simulate = true;
    simulationSettings.plot = true;
    simulationSettings.animate = false;
        simulationSettings.pauseStart = false;
    simulationSettings.stepTime = 0.5;              % s
    simulationSettings.T = 5;                       % s
    simulationSettings.dt = 0.01;                   % s 
    
% model settings
simulationSettings.measurementHeight = 0.95;    % []*L
simulationSettings.dampingModel = 'modal';      % 'Modal','proportional' or 'none'
simulationSettings.Input = 'force';             % 'force' for force input, or 'disp' for disp input   
    simulationSettings.Kbase = 1e10;            % only in 'disp' mode
simulationSettings.Nsensors = 4;                % number of sensor patches
simulationSettings.Nactuators = 4;              % number of actuators
        
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
plotSettings.links = false;                     % Plot the links (edges)
    plotSettings.linkNumbers = false;           % link numbers
plotSettings.sensor = true;                     % plot the sensor
plotSettings.sensorPlot = true;                 % plot the sensor output
plotSettings.base = true;                       % Plot the location of the base in the measurement


%% Build FEM model of the beam
disp('Running!-------------------------------------------------------------')
[FEM,beam,PDEbeam] = buildModel(beam,sensor,actuator,plotSettings,simulationSettings);
numNodes = length(beam.nodes);
numLinks = length(beam.links);
numFaces = length(beam.faces);

K = FEM.K;
M = FEM.M;

% Get modal results
wmin = -0.1;    % -0.1 rad/s to capture everything from 0 ->

% Modal decomposition! P395
modal = solve(PDEbeam,'FrequencyRange',[wmin,wmax]);   % Solve
Nmodes = length(modal.NaturalFrequencies);

% Find node to sense!
nodesxy = PDEbeam.Mesh.Nodes;
measurey = simulationSettings.measurementHeight*beam.L;
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
Bd = zeros(numNodes*2,2);                                       % External Input matrix in d space
Bd(1,2) = 1;                                                    % Force on node 1 in x; (Force input)                                                
Bd(2,2) = 1;                                                    % Force on node 2 in x; (Force input)
Bd(3,1) = 1;                                                    % Force on node 3 in x; (Base Force)
Bd(4,1) = 1;                                                    % Force on node 4 in x, (Base Force)

% Output
Cd = zeros(2,numNodes*2);                                       % Measurement matrix in d space
Cd(1,interpNodes(:,1)) = [1-alpha;alpha];                       % Measure at selected point
Cd(2,1) = 1;                                                    % Measure base displacement


%% modal system modelling
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
        SS_modal = SS_modal(1,[1,2]);
        
    case {'disp'}
        % To get the input of the system bo be the displacement of the base
        % of the beam, it needs to be fixed by a stiff spring. This is done
        % by applying state feedback with a gain representing a spring
        % constant. 

        K_spring = ss(simulationSettings.Kbase);    % get stiffness matrix

        SS_modal = feedback(SS_modal,K_spring,2,2,-1);
        DCgain = dcgain(SS_modal(1,2));
        SS_modal.B = [B_modal(:,1),B_modal(:,2)/DCgain];

        SS_modal = SS_modal(1,[1,2]);               % From [F,r] to sensor out. 
        SS_modal.Inputname = {'Force in','base Ref'};
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

%% Time simulation
if simulationSettings.simulate ==  true
    disp('Simulating...')
    tvec = 0:simulationSettings.dt:simulationSettings.T;
    
    % Time simulation!
    q0 = zeros(Nmodes*2,1);
    Asim = SS_modal.A;
    Bsim = SS_modal.B;
    [t,y] = ode45(@(t,z) sysFun(t,z,Asim,Bsim,simulationSettings),tvec,q0);
    z = y';
    output = SS_modal.C*z;
    baseLocation = [Cphi(2,:),zeros(1,Nmodes)]*z;
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
        measurement = zeros(2,length(t));        
        
        if simulationSettings.animate == true
            Nsteps = length(t);
        else
            Nsteps = 1;
        end

        % Run animation loop! (Only once if animate is turned off)
        for i = 1:Nsteps
            tstart = tic;
            d = Phi*z(1:Nmodes,i);

            % Update Beam
            beam.update(d);
            
            % Delete old plot
            if ~isempty(simPlots)
                delete(simPlots)
            end

            % Plot beam
            simPlots = plotModel(beam,plotSettings,simulationSettings,beamAx); 
            timeText = text(beamAx, 0,beam.L*1.1,['Time: ',num2str(round(t(i),1)),'/',num2str(t(end)),' s'],...
                                                        'HorizontalAlignment','center');
            simPlots = [simPlots, timeText];
            
            % Plot laser
            if plotSettings.sensor == true
                laser = plot(beamAx,[beamAx.XLim(1) output(i)-beam.h/2],[1,1]*simulationSettings.measurementHeight*beam.L,'r','lineWidth',2);
                simPlots = [simPlots, laser];
            end
            
            % Plot laser measurement
            if plotSettings.sensorPlot == true
                if simulationSettings.animate == true
                    if i > 1
                        xs = [t(i-1),t(i)];
                        ys = [output(i-1),output(i)];
                        if plotSettings.base == true
                            yb = [baseLocation(i-1),baseLocation(i)];
                        end
                    else
                        xs = t(i);
                        ys = output(i);
                        yb = baseLocation(i);
                    end
                else
                    xs = t;
                    ys = output;
                    yb = baseLocation;
                end
                    if plotSettings.base == true
                        basePlot = plot(measurementAx,xs,yb,'r');
                    end
                    sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
            end
            
            % Plot input and output as dots
            inputNodes = find(Bd);
            inputNode = inputNodes(1);
            
            inputN = plot(beamAx,beam.nodes{inputNode}.pos(1),beam.nodes{inputNode}.pos(2),'g.','MarkerSize',12);
            simPlots = [simPlots, inputN];
            
            alpha = simulationSettings.interpPoint(1,3);
            interpNodes = simulationSettings.interpPoint(:,2);
            
            outputNode = interpNodes(1)+alpha*(interpNodes(2)-interpNodes(1));
            outputN = plot(beamAx,output(i)-beam.h/2,outputNode,'r.','MarkerSize',12);
            simPlots = [simPlots, outputN];
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
            updateBeam(links,nodes,faces,A*(Phi(:,k)/norm(Phi(:,k))));
            % Plot Beam
            plotBeam(nodes,links,faces,plotSettings,simulationSettings,ax);
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
            else
                F = 0;
            end
        case {'disp'}
            F = 0;
            if t > simulationSettings.stepTime
                r = 1;
            else
                r = 0;
            end
    end        
    % Evaluate system
    qd = A*q+B*[F,r]';
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
              
              
              