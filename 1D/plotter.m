function simPlots = plotter(model,m)

% Check if model is simulated
if isempty(model.simulationData)
    warning("Model has not been simulated yet, can't plot response")
    return
end

% Set alpha based on model
alphas = linspace(model.plotSettings.alphaMin,1,model.simulationSettings.nModels);

modelNumber = model.number;
Alpha = alphas(modelNumber);

% Defenitions
L =         model.modelSettings.L;          

simulationSettings =    model.simulationSettings;
modelSettings =         model.modelSettings;
plotSettings =          model.plotSettings;

qfull =     model.simulationData.qfull(:,:,m);             % full states in q space
y =         model.simulationData.yfull(:,:,m);
Udist =     model.simulationData.Udist(:,:,m);

y_MF = model.simulationData.yfull_MF(:,:,m);
y_LO = model.simulationData.yfull_LO(:,:,m);
y_KF = model.simulationData.yfull_KF(:,:,m);
y_AKF = model.simulationData.yfull_AKF(:,:,m);
y_DKF = model.simulationData.yfull_DKF(:,:,m);
y_GDF = model.simulationData.yfull_GDF(:,:,m);

u_DKF = model.simulationData.ufull_DKF(:,:,m);
u_GDF = model.simulationData.ufull_GDF(:,:,m);

nPatches = length(model.modelSettings.patches);
nAcc = length(modelSettings.Acc);
Nmodes = modelSettings.Nmodes;
t = 0:simulationSettings.dt:simulationSettings.T;

% Setting up figure
if isempty(findobj('Name','Simulation results'))
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
        beamAx.Tag = 'beamAx';
    subplot(12,3,[2,3,5,6])
        hold on
        grid on
        ylabel 'displasement [m]'
        title 'Laser Measurement'
        measurementAx = gca;
        measurementAx.Tag = 'measurementAx';
        xlim([0,simulationSettings.T]);
        if min(y(1,:)) ~= max(y(1,:))
            ylim(1.1*[min(y(1,:)),max(y(1,:))])
        end
    subplot(12,3,[11,12,14,15]);
        hold on
        grid on
        switch model.modelSettings.gauges
            case {'piezo'}
                ylabel 'V'
                title 'Piezo measurement' 
            case {'resistive'}
                ylabel '$$\epsilon$$'
                title 'Strain measurement' 
        end
    
        piezoAx = gca;
        piezoAx.Tag = 'piezoAx';
        xlim([0,simulationSettings.T]);
    subplot(12,3,[20,21,23,24])
        hold on
        grid on
        xlabel 'Time [s]'
        ylabel '$m/s^2$'
        title 'Accelometer Outputs'   
        accAx = gca;
        accAx.Tag = 'accAx';
        xlim([0,simulationSettings.T]);
    subplot(12,3,28:36)
        bodeAx = gca;
        bodeAx.Tag = 'bodeAx';
    
    movegui(a,'northwest')
else
    % Find figures again. 
    beamAx = findobj('Tag','beamAx');
    measurementAx = findobj('Tag','measurementAx');
    piezoAx = findobj('Tag','piezoAx');
    accAx = findobj('Tag','accAx');
    bodeAx = findobj('Tag','bodeAx');
end
    
%% Bode plot!
if modelNumber == 1
    [~] = model.showBode(bodeAx);
end

%% Beam plot!
i = 1;  % For now only plot initial orientation possibility to animate still there. 
q = model.simulationData.qfull(:,i); % State to be plotted
simPlots = model.showBeam(beamAx,q);

%% Output measurement plot!
plot(measurementAx,t,y(1,:),'k'); % Sensor
if any(ismember(simulationSettings.observer,'MF'))
    plot(measurementAx,t,y_MF(1,:),'color',[0 0.4470 0.7410 Alpha]);
end
if any(ismember(simulationSettings.observer,'LO'))
    plot(measurementAx,t,y_LO(1,:),'color',[0.8500 0.3250 0.0980 Alpha]);
end
if any(ismember(simulationSettings.observer,'KF'))
    plot(measurementAx,t,y_KF(1,:),'color',[0.9290 0.6940 0.1250 Alpha]);
end
if any(ismember(simulationSettings.observer,'AKF'))
    plot(measurementAx,t,y_AKF(1,:),'color',[0.4940 0.1840 0.5560 Alpha]);
end
if any(ismember(simulationSettings.observer,'DKF'))
    plot(measurementAx,t,y_DKF(1,:),'color',[0.3010 0.7450 0.9330 Alpha]);
end
if any(ismember(simulationSettings.observer,'GDF'))
    plot(measurementAx,t,y_GDF(1,:),'color',[0.6350 0.0780 0.1840 Alpha]);
end

legend(measurementAx,['True' simulationSettings.observer])

%% Piezo plot
colors = [1 0 0 Alpha; 0 1 0 Alpha; 0 0 1 Alpha; 0 0 0 Alpha];
lineStyles = ['-','--'];

if nPatches > 0                    % If there are patches
    piezoOutputs = y(2:nPatches+1,:);
    for i = 1:nPatches
        if i > length(colors)
            j = i-length(colors);
        else
            j = i;
        end

        color = colors(j,:);
        lineStyle = lineStyles(ceil(i/length(colors)));
        plot(piezoAx,t,piezoOutputs(i,:),'Color',color,'LineStyle',lineStyle);
    end
end

%% Accelerometer plot
if nAcc > 0                         % If there are accelerometers
    ts = ones(nAcc,1)*t;
    accOutputs = y(nPatches+2:end,:);
    for i = 1:nAcc
        if i > length(colors)
            j = i-length(colors);
        else
            j = i;
        end

        color = colors(j,:)';
        lineStyle = lineStyles(ceil(i/length(colors)));
        plot(accAx,ts',accOutputs(i,:)','Color',color,'LineStyle',lineStyle);
    end
end



%% State plot
if plotSettings.statePlot == true
    PNmodes = plotSettings.states; % If you only want to plot the first Nmodes modes. 
    axes = gobjects(PNmodes);
    
    name = 'Observer and true states';
    if isempty(findobj('Name',name))
        b = figure('Name',name);
        sgtitle('Modal states')
        hold on
        for i = 1:PNmodes
            subname = ['State ', i];

            subplot(PNmodes,1,i)
            hold on
            grid on
            ylabel(['$\eta$',num2str(i)])  
            xlim([0,simulationSettings.T]);
            if min(qfull(i,:)) ~= max(qfull(i,:))
                ylim(1.1*[min(qfull(i,:)),max(qfull(i,:))])  
            else
            end
            axes(i) = gca;
            axes(i).Tag = subname;
            plot(gca,t,qfull(i,:),'k');             % Plot true states
    
            movegui(b,'north')
        end
    else
        for i = 1:PNmodes
            subname = ['State ', i];
            axes(i) = findobj('Tag',subname);
        end
    end
    
    % Plot modal filter
    if any(ismember(simulationSettings.observer,'MF'))
        qfull_MF = model.simulationData.qfull_MF;
        for i = 1:PNmodes
            plot(axes(i),t,qfull_MF(i,:,m),'color',[0 0.4470 0.7410 Alpha]);
        end
    end

    % Plot Luenberger observer
    if any(ismember(simulationSettings.observer,'LO'))
        qfull_LO = model.simulationData.qfull_LO;
        for i = 1:PNmodes
            plot(axes(i),t,qfull_LO(i,:,m),'color',[0.8500 0.3250 0.0980 Alpha]);
        end
    end 
    
    % Plot conventional Kalman Filter
    if any(ismember(simulationSettings.observer,'KF'))
        qfull_KF = model.simulationData.qfull_KF;
        for i = 1:PNmodes
            plot(axes(i),t,qfull_KF(i,:,m),'color',[0.9290 0.6940 0.1250 Alpha]);
        end
    end

    % Plot Augmented Kalman Filter
    if any(ismember(simulationSettings.observer,'AKF'))
        qfull_AKF = model.simulationData.qfull_AKF;
        for i = 1:PNmodes
            plot(axes(i),t,qfull_AKF(i,:,m),'color',[0.4940 0.1840 0.5560 Alpha]);
        end
    end

    % Plot Dual Kalman Filter
    if any(ismember(simulationSettings.observer,'DKF'))
        qfull_DKF = model.simulationData.qfull_DKF;
        for i = 1:PNmodes
            plot(axes(i),t,qfull_DKF(i,:,m),'color',[0.3010 0.7450 0.9330 Alpha]);
        end

    end

    % Plot Giljins de Moor Filter (GDF)
    if any(ismember(simulationSettings.observer,'GDF'))
        qfull_GDF = model.simulationData.qfull_GDF;
        for i = 1:PNmodes
            plot(axes(i),t,qfull_GDF(i,:,m),'color',[0.6350 0.0780 0.1840 Alpha]);
        end
    end
    legend(axes(1),['True' simulationSettings.observer])
end


%% Input plot
if plotSettings.inputSequence == true
    if isempty(findobj('Name','System input sequence'))
        c = figure('Name' ,'System input sequence');
        hold on
        grid on
        ylabel 'Force input [N]'
        xlabel 'Time [s]'
        title 'System input sequence'
        inputAx = gca;
        inputAx.Tag = 'inputAx';
        xlim([0,simulationSettings.T]);
        plot(inputAx,t,Udist,'k'); 
    
        movegui(c,'northeast');
    else
        inputAx = findobj('Tag','inputAx');
    end

    % Plot Augmented Kalman Filter
    if any(ismember(simulationSettings.observer,'AKF'))
        plot(inputAx,t,qfull_AKF(Nmodes*2+1,:,m),'color',[0.4940 0.1840 0.5560 Alpha]);
    end
    
    % Plot Dual Kalman Filter
    if any(ismember(simulationSettings.observer,'DKF'))
        plot(inputAx,t,u_DKF(1,:,m),'color',[0.3010 0.7450 0.9330 Alpha]);
    end

    % Plot Giljins de Moor Filter (GDF)
    if any(ismember(simulationSettings.observer,'GDF'))
        plot(inputAx,t,u_GDF(1,:,m),'color',[0.6350 0.0780 0.1840 Alpha]);
    end
    
    AKFLocation = ismember(simulationSettings.observer,"AKF");
    DKFLocation = ismember(simulationSettings.observer,"DKF");
    GDFLocation = ismember(simulationSettings.observer,"GDF");
    theRest = logical(AKFLocation + DKFLocation + GDFLocation);
    legend(inputAx,['True' simulationSettings.observer(theRest)])
end

%% Mode shape plot
if plotSettings.modes > 0
    d = figure('Name','Mode Shapes');
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