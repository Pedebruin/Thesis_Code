function simPlots = plotter(model)
disp('Plotting simulation...')

% Check if model is simulated
if isempty(model.simulationData)
    warning("Model has not been simulated yet, can't plot response")
    return
end

% Unpacking and defining for readability
L =         model.modelSettings.L;  
numEl =     model.numEl;
numNodes =  model.numNodes;
Nmodes =    model.modelSettings.Nmodes;
Phi =       model.Phi;
elements =  model.elements;
sys =       model.sys;          

simulationSettings =    model.simulationSettings;
modelSettings =         model.modelSettings;
plotSettings =          model.plotSettings;

qfull =     model.simulationData.qfull;             % full states in q space
y =         model.simulationData.y;
Udist =     model.simulationData.Udist;

nPatches = length(modelSettings.sElements)/modelSettings.nsElementsP;
nAcc = length(modelSettings.Acc);
t = 0:simulationSettings.dt:simulationSettings.T;
dfull = [Phi, zeros(numNodes*2,Nmodes);         % full states in d space
        zeros(numNodes*2,Nmodes),Phi]*qfull;

% Setting up figure
figure('Name','Simulation results')
subplot(12,3,(0:8)*3+1) % Beam plot
    hold on
    grid on
    xlabel 'm'
    ylabel 'm'
    axis equal
    xlim([-L/6,L/6]);
    ylim([0,1.2*L])
    title 'Beam plot'
    beamAx = gca;
subplot(12,3,[2,3,5,6])
    hold on
    grid on
    ylabel 'displasement [m]'
    title 'Laser Measurement'
    measurementAx = gca;
    xlim([0,simulationSettings.T]);
subplot(12,3,[11,12,14,15]);
    hold on
    grid on
    ylabel 'V'
    title 'Piezo outputs'   
    piezoAx = gca;
    xlim([0,simulationSettings.T]);
subplot(12,3,[20,21,23,24])
    hold on
    grid on
    xlabel 'time [s]'
    ylabel '$m/s^2$'
    title 'Accelometer Outputs'   
    accAx = gca;
    xlim([0,simulationSettings.T]);
subplot(12,3,28:36)
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
plotoptions.grid = 'on';
plotoptions.PhaseWrapping = 'off';

bodeplot(bodeAx,sys(1,1),plotoptions);

% Plot analytical omega1
xline(bodeAx,model.analyticalf1)

%% Beam plot!
%%%------------------------------------------------------------------------
simPlots = [];     
i = 1;  % For now only plot initial orientation possibility to animate still there. 

% update elements
for j = 1:length(elements)
    elements(j).update(dfull(:,i));  
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

timeText = text(beamAx, 0,L*1.1,['T= ',num2str(round(t(i),1)),'/',num2str(t(end)),' s'],...
                                            'HorizontalAlignment','center');
simPlots = [simPlots, timeText];

%%%------------------------------------------------------------------------
% Plot laser
if plotSettings.sensor == true
    laserx = y(1,i);
    laser = plot(beamAx,[beamAx.XLim(1) laserx],[1,1]*modelSettings.measurementHeight*L,'r','lineWidth',2);
    simPlots = [simPlots, laser];
end

%% Laser plot
plot(measurementAx,t,y(1,:),'color',[0.3010 0.7450 0.9330]); % Sensor

if plotSettings.Input == true
    plot(measurementAx,t,Udist,'r'); % Sensor
end

%% Piezo plot
if nPatches > 0                    % If there are patches
    piezoOutputs = y(2:nPatches+1,:);
    plot(piezoAx,t,piezoOutputs);
end

%% Accelerometer plot
if nAcc > 0                         % If there are accelerometers
    ts = ones(nAcc,1)*t;
    accOutputs = y(nPatches+2:end,:);
    plot(accAx,ts',accOutputs');
end

drawnow;

if plotSettings.statePlot == true
    Nmodes = Nmodes; % If you only want to plot the first Nmodes modes. 
    if any(ismember(simulationSettings.observer,'LO'))
        qfull_LO = model.simulationData.qfull_LO;
        figure('Name','LO')
        sgtitle('Leuenberger Observer')
        
        for i = 1:Nmodes
            subplot(Nmodes,1,i)
            hold on
            grid on
            ylabel(['$\eta$',num2str(i)])  
            xlim([0,simulationSettings.T]);

            plot(gca,t,qfull(i,:));
            plot(gca,t,qfull_LO(i,:));
        end
    end 

    if any(ismember(simulationSettings.observer,'AKF'))
        qfull_AKF = model.simulationData.qfull_AKF;
        figure('Name','AKF')
        sgtitle('Augmented Kalman Filter')
        for i = 1:Nmodes
            subplot(Nmodes,1,i)
            hold on
            grid on
            ylabel(['$\eta$',num2str(i)])  
            xlim([0,simulationSettings.T]);
    
            plot(gca,t,qfull(i,:));
            plot(gca,t,qfull_AKF(i,:));
        end
    end
    if any(ismember(simulationSettings.observer,'DKF'))
        qfull_DKF = model.simulationData.qfull_DKF;
        figure('Name','DKF')
        sgtitle('Dual Kalman Filter')
        for i = 1:Nmodes
            subplot(Nmodes,1,i)
            hold on
            grid on
            ylabel(['$\eta$',num2str(i)])  
            xlim([0,simulationSettings.T]);
    
            plot(gca,t,qfull(i,:));
            plot(gca,t,qfull_DKF(i,:));
        end

    end
    if any(ismember(simulationSettings.observer,'GDF'))
        qfull_GDF = model.simulationData.qfull_GDF;
        figure('Name','GDF')
        sgtitle('Giljins, de Moor Filter')
        for i = 1:Nmodes
            subplot(Nmodes,1,i)
            hold on
            grid on
            ylabel(['$\eta$',num2str(i)])  
            xlim([0,simulationSettings.T]);
    
            plot(gca,t,qfull(i,:));
            plot(gca,t,qfull_GDF(i,:));
        end
    end
    drawnow;

end
end