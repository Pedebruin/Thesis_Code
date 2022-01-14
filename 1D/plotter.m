function simPlots = plotter(model)

% Check if model is simulated
if isempty(model.simulationData)
    warning("Model has not been simulated yet, can't plot response")
    return
end
L =         model.modelSettings.L;          

simulationSettings =    model.simulationSettings;
modelSettings =         model.modelSettings;
plotSettings =          model.plotSettings;

qfull =     model.simulationData.qfull;             % full states in q space
y =         model.simulationData.yfull;
Udist =     model.simulationData.Udist;

y_MF = model.simulationData.yfull_MF;
y_LO = model.simulationData.yfull_LO;
y_KF = model.simulationData.yfull_KF;
y_AKF = model.simulationData.yfull_AKF;
y_DKF = model.simulationData.yfull_DKF;
y_GDF = model.simulationData.yfull_GDF;

u_DKF = model.simulationData.ufull_DKF;
u_GDF = model.simulationData.ufull_GDF;

nPatches = length(model.modelSettings.patches);
nAcc = length(modelSettings.Acc);
t = 0:simulationSettings.dt:simulationSettings.T;

% Setting up figure
a = figure('Name','Simulation results');
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
    if min(y(1,:)) ~= max(y(1,:))
        ylim(1.1*[min(y(1,:)),max(y(1,:))])
    else
    end
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

movegui(a,'northwest')
    
%% Bode plot!
[~] = model.showBode(bodeAx);

%% Beam plot!
i = 1;  % For now only plot initial orientation possibility to animate still there. 
         
q = model.simulationData.qfull(:,i); % State to be plotted

simPlots = model.showBeam(beamAx,q);

%% Output measurement plot!
plot(measurementAx,t,y(1,:),'k'); % Sensor
if any(ismember(simulationSettings.observer,'MF'))
    plot(measurementAx,t,y_MF(1,:),'color',[0 0.4470 0.7410]);
end
if any(ismember(simulationSettings.observer,'LO'))
    plot(measurementAx,t,y_LO(1,:),'color',[0.8500 0.3250 0.0980]);
end
if any(ismember(simulationSettings.observer,'KF'))
    plot(measurementAx,t,y_KF(1,:),'color',[0.9290 0.6940 0.1250]);
end
if any(ismember(simulationSettings.observer,'AKF'))
    plot(measurementAx,t,y_AKF(1,:),'color',[0.4940 0.1840 0.5560]);
end
if any(ismember(simulationSettings.observer,'DKF'))
    plot(measurementAx,t,y_DKF(1,:),'color',[0.3010 0.7450 0.9330]);
end
if any(ismember(simulationSettings.observer,'GDF'))
    plot(measurementAx,t,y_GDF(1,:),'color',[0.6350 0.0780 0.1840]);
end

legend(measurementAx,['True' simulationSettings.observer])

%% Piezo plot
if nPatches > 0                    % If there are patches
    piezoOutputs = y(2:nPatches+1,:);
    plot(piezoAx,t,piezoOutputs);
end

piezoAx.ColorOrder = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
piezoAx.LineStyleOrder = {'-','--'};


%% Accelerometer plot
if nAcc > 0                         % If there are accelerometers
    ts = ones(nAcc,1)*t;
    accOutputs = y(nPatches+2:end,:);
    plot(accAx,ts',accOutputs');
end

accAx.ColorOrder = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
accAx.LineStyleOrder = {'-','--'};


%% State plot
if plotSettings.statePlot == true
    Nmodes = plotSettings.states; % If you only want to plot the first Nmodes modes. 
    axes = gobjects(Nmodes);
    
    b = figure('Name','Observer and true states');
    sgtitle('Modal states')
    hold on
    for i = 1:Nmodes
        subplot(Nmodes,1,i)
        hold on
        grid on
        ylabel(['$\eta$',num2str(i)])  
        xlim([0,simulationSettings.T]);
        if min(qfull(i,:)) ~= max(qfull(i,:))
            ylim(1.1*[min(qfull(i,:)),max(qfull(i,:))])  
        else
        end
        axes(i) = gca;
        plot(gca,t,qfull(i,:),'k');             % Plot true states

        movegui(b,'north')
    end
    
    % Plot modal filter
    if any(ismember(simulationSettings.observer,'MF'))
        qfull_MF = model.simulationData.qfull_MF;
        for i = 1:Nmodes
            plot(axes(i),t,qfull_MF(i,:),'color',[0 0.4470 0.7410]);
        end
    end

    % Plot Luenberger observer
    if any(ismember(simulationSettings.observer,'LO'))
        qfull_LO = model.simulationData.qfull_LO;
        for i = 1:Nmodes
            plot(axes(i),t,qfull_LO(i,:),'color',[0.8500 0.3250 0.0980]);
        end
    end 
    
    % Plot conventional Kalman Filter
    if any(ismember(simulationSettings.observer,'KF'))
        qfull_KF = model.simulationData.qfull_KF;
        for i = 1:Nmodes
            plot(axes(i),t,qfull_KF(i,:),'color',[0.9290 0.6940 0.1250]);
        end
    end

    % Plot Augmented Kalman Filter
    if any(ismember(simulationSettings.observer,'AKF'))
        qfull_AKF = model.simulationData.qfull_AKF;
        for i = 1:Nmodes
            plot(axes(i),t,qfull_AKF(i,:),'color',[0.4940 0.1840 0.5560]);
        end
    end

    % Plot Dual Kalman Filter
    if any(ismember(simulationSettings.observer,'DKF'))
        qfull_DKF = model.simulationData.qfull_DKF;
        for i = 1:Nmodes
            plot(axes(i),t,qfull_DKF(i,:),'color',[0.3010 0.7450 0.9330]);
        end

    end

    % Plot Giljins de Moor Filter (GDF)
    if any(ismember(simulationSettings.observer,'GDF'))
        qfull_GDF = model.simulationData.qfull_GDF;
        for i = 1:Nmodes
            plot(axes(i),t,qfull_GDF(i,:),'color',[0.6350 0.0780 0.1840]);
        end
    end
    legend(axes(1),['True' simulationSettings.observer])
end


%% Input plot
if plotSettings.inputSequence == true
    c = figure('Name' ,'System input sequence');
    hold on
    grid on
    ylabel 'force Input [N]'
    title 'System input sequence'
    inputAx = gca;
    xlim([0,simulationSettings.T]);
    plot(inputAx,t,Udist,'k'); 

    movegui(c,'northeast');

    % Plot Augmented Kalman Filter
    if any(ismember(simulationSettings.observer,'AKF'))
        plot(inputAx,t,qfull_AKF(Nmodes*2+1,:),'color',[0.4940 0.1840 0.5560]);
    end
    
    % Plot Dual Kalman Filter
    if any(ismember(simulationSettings.observer,'DKF'))
        plot(inputAx,t,u_DKF(1,:),'color',[0.3010 0.7450 0.9330]);
    end

    % Plot Giljins de Moor Filter (GDF)
    if any(ismember(simulationSettings.observer,'GDF'))
        plot(inputAx,t,u_GDF(1,:),'color',[0.6350 0.0780 0.1840]);
    end
    
    AKFLocation = ismember(simulationSettings.observer,"AKF");
    DKFLocation = ismember(simulationSettings.observer,"DKF");
    GDFLocation = ismember(simulationSettings.observer,"GDF");
    theRest = logical(AKFLocation + DKFLocation + GDFLocation);
    legend(inputAx,['True' simulationSettings.observer(theRest)])
end

%% Mode shape plot
if plotSettings.modes > 0
    d = figure('Name',['Mode Shapes']);
    sgtitle('Mode shapes')
    for i = 1:plotSettings.modes
        subplot(1,plotSettings.modes,i)
            hold on
            grid on
            xlabel 'm'
            ylabel 'm'
            axis equal
            xlim([-L/6,L/6]);
            ylim([0,1.2*L])
            title(['Mode ',num2str(i)])
            simPlots = model.showMode(gca,i);
    end
    movegui(d,"southwest");
end
drawnow;
end