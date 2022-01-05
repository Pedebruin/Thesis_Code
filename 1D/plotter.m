function plots = plotter(y,y_fb,dmat,sys,sys_fb,elements,analyticalf1,plotSettings,modelSettings,simulationSettings)
L = modelSettings.L;  
numEl = length(elements);

disp('Plotting simulation...')
    figure()
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
    subplot(12,3,[2,3,5,6,8,9])
        hold on
        grid on
        ylabel 'displasement [m]'
        title 'Laser Measurement'
        measurementAx = gca;
        xlim([0,simulationSettings.T]);
    subplot(12,3,[11,12,14,15,17,18]);
        hold on
        grid on
        xlabel 'time [s]'
        ylabel 'V'
        title 'Output voltages'   
        stateAx = gca;
        xlim([0,simulationSettings.T]);
    subplot(12,3,[20,21,23,24,26,27])
        hold on
        grid on
        xlabel 'time [s]'
        ylabel '$m/s^2$'
        title 'Output accelerations'   
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
    plotOptions.grid = 'on';
    plotOptions.PhaseWrapping = 'off';
    plotOptions.XLimMode = {'manual','manual'};

    if modelSettings.fb == true
        bod = bodeplot(bodeAx,sys(1,1),sys_fb(1,1),plotoptions);
    else
        bod = bodeplot(bodeAx,sys(1,1),plotoptions);
    end
    
    % Plot analytical omega1
    xline(bodeAx,analyticalf1)
  
    simPlots = [];     
    
    Nsteps = 1;

    % Run animation loop! (Only once if animate is turned off)
    for i = 1:Nsteps
        tstart = tic;

        d = dmat(:,end);


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
        
        % Plot laser measurement, voltages and accelerations
        if plotSettings.sensorPlot == true
            xs = t;
            baseLocation = dmat(1,:);
            ys = measurement;
            if nsElements > 0                   % If there are patches
                ysbeam = y(3:nPatches+2,:);
                if modelSettings.fb == true      
                    ys_fb = measurement_fb;
                end
            end
            sensorPlot = plot(measurementAx,xs,ys,'color',[0.3010 0.7450 0.9330]);
            if nsElements > 0                   % If there are patches
                statePlot = plot(stateAx,xs,ysbeam);
                if modelSettings.fb == true
                    sensorfbPlot = plot(measurementAx,xs,ys_fb,'color',[0.8500 0.3250 0.0980]);
                end
            end
        end
        
        drawnow;
        
        % Measure elapsed time and delay
        elapsed = toc(tstart);
        
        if simulationSettings.animate == true
            %pause(max(simulationSettings.dt-elapsed,0));
        end
        if simulationSettings.pauseStart == true && i == 1 && simulationSettings.animate == true
            disp('PAUSED, Press any key to continue');
            pause;
        end
    end
end