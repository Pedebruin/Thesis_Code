function plotDataset(dataset)
    close all
    addpath('./dataSets')
    set(0,'defaultTextInterpreter','latex');
    load(dataset)
    
    t = data{5}.Values.Time;
    
    [Beam, sBeam, modelSettings,...
    simulationSettings, plotSettings,simulinkSettings...
    ,LO,KF,AKF,DKF,GDF,~, ~] = configure(0);
    
    positionMeasurement = data{5}.Values.Data';
    reference = data{3}.Values.Data';
    Effort = data{4}.Values.Data';
    strainMeasurement = data{2}.Values.Data'*modelSettings.strainCorrectionGain;
    strainMeasurement = strainMeasurement - mean(strainMeasurement); % Because I didn't want to calibrate every time
    
    accelerationMeasurement = double(data{1}.Values.Data')/simulationSettings.accSensitivity;
    accelerationMeasurement = accelerationMeasurement-mean(accelerationMeasurement); % Needs to be solved for real-time application!
    accelerationMeasurement = lowpass(accelerationMeasurement,1e3,1/simulationSettings.dt);
    
    figure()
    hold on
    subplot(4,5,[1:5])
        hold on
%         title('Position and reference')
        title('Position')
        ylabel('[m]')
        plot(t,reference,'Color',[0.8500 0.3250 0.0980])
        plot(t,positionMeasurement,'Color',[0 0.4470 0.7410])
        legend({'Reference','Position'})
    subplot(4,5,[6:10])
        hold on
        title('Acceleration')
        ylabel('[m/s$^2$]')
        plot(t,accelerationMeasurement)
    subplot(4,5,[11:15])
        hold on
        title('Strain')
        ylabel('[m/m]')
        plot(t,strainMeasurement)
    subplot(4,5,[16:20])
        hold on
        title('Control effort')
        xlabel('time [s]')
        ylabel('[$\%$]')
        plot(t,Effort/100)
end