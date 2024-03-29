clear
close all
stored = fileparts(which('sysID.m'));
chdir(stored)
addpath('./../c2000/Experiments/')
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15,'defaultLegendFontSize',15);  

% Bode settings
P = bodeoptions;
P.FreqUnits = 'Hz';
P.PhaseMatching = 'on';
P.PhaseMatchingValue = 180;
P.PhaseMatchingFreq = 260*2*pi;
P.Grid = 'on';

%% See Influence of geometric parameters! (Optional offcourse)
% Plot parameter variations!
varyParameters = false;
if varyParameters
    N = 10;                                 % Amount of bode plots?
    evaluateParameters([],"L",[85e-3,110e-3],N,"Ssteel",2,1);
    evaluateParameters([],"b",[10e-3,20e-3],N,"Ssteel",2,1);
    evaluateParameters([],"h",[0.1e-3,1e-3],N,"Ssteel",2,1);
%     evaluateParameters([],"E",[68.5E9,210E9],N,"Ssteel",2,1);
%     evaluateParameters([],"rho",[2700,7800],N,"Ssteel",2,1);
    evaluateParameters([],"stageMass",[350e-3/4,650e-3/4],N,"Ssteel",2,1);
end

% Look at specific cases
specificConfigurations = false;
if specificConfigurations
    evaluateParameters([],[],[],1,"Aluminium",1,1); % Aluminium standard beam
    evaluateParameters(gca,["h"],[0.1e-3,1e-3],10,"Ssteel",1,1); % Spring steel ISI 1095 with width variation
    evaluateParameters(gca,["h"],0.4e-3,1,"Ssteel",1,1); % Specific width
end

%% Optimise patch length and placement!
optimisePatches = false;
if optimisePatches
    evaluateParameters([],"patchL",[15e-3,30e-3],10,"Ssteel",2,1); % Show all 
    X1 = fminsearch(@optfnc,[75e-3,22]*1e-3);
    X2 = fminsearch(@optfnc,[1,80]*1e-3);
    
%     evaluateParameters(gca,["patches","patchL"],X,1,"Ssteel",2,1); % Plot optimal config
%     evaluateParameters([],["patches","patchL",],X,1,"Ssteel",1,1);
    
    gridSize = 50;
    L = 100e-3;
    lengths = linspace(1e-3,99e-3,gridSize)';
    locations = linspace(0,0.99,gridSize);
    
    grid = zeros(gridSize,gridSize);
    for i = 1:gridSize
        for j = 1:gridSize
            grid(i,j) = optfnc([locations(i),lengths(j)]);
        end
    end
    
    figure()
    hold on
        % Plot response function
    subplot(4,4,[1:3,5:7,9:11,13:15])
        hold on
        title('Second mode observability as a function of patch length and position')
        ylabel('Patch Position [m]')
        xlabel('Patch Length [m]')
        zlabel('Peak Magnitude []')
        surface = gca;
        surf(surface,ones(gridSize,1)*locations*L,lengths*ones(1,gridSize),-grid)
    
        plot3(surface,X1(2),X1(1)*L,-optfnc(X1),'ok','MarkerFaceColor','k')
        text(surface,X1(2),X1(1)*L,-optfnc(X1)+2e-7,['1: [',num2str(round(X1(2),3)),',',num2str(round(X1(1)*L,3)),']'],...
            'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        plot3(surface,X2(2),X2(1)*L,-optfnc(X2),'ok','MarkerFaceColor','k')
        text(surface,X2(2),X2(1)*L,-optfnc(X2)+2e-7,['2: [',num2str(round(X2(2),3)),',',num2str(round(X2(1)*L,3)),']'],...
            'HorizontalAlignment','center','VerticalAlignment','bottom')

        % Plot 1st and 2nd optimum
     subplot(4,4,[4,8])   
        hold on
        title('Optimum 1')
        opt1 = gca;
        [Beam, sBeam, modelSettings,...
                simulationSettings, plotSettings,~,~,~,~,~,~,~,~] = configure(1,["patches","patchL"],X1,"Ssteel");
        
        plotSettings.inputForce = false;
        plotSettings.sensor = false;

        SYS = model(); % create model object
        [SYS,~] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 
        SYS.showBeam(opt1,zeros(SYS.nq,1));

    subplot(4,4,[12,16])
        hold on
        title('Optimum 2')
        opt2 = gca;
        [Beam, sBeam, modelSettings,...
                simulationSettings, plotSettings,~,~,~,~,~,~,~,~] = configure(1,["patches","patchL"],X2,"Ssteel");
        
        plotSettings.inputForce = false;
        plotSettings.sensor = false;

        SYS = model(); % create model object
        [SYS,~] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 
        SYS.showBeam(opt2,zeros(SYS.nq,1));   
end
    
%% Configure using standard settings
identify = true;
if identify == true
    [Beam, sBeam, modelSettings,...
            simulationSettings, plotSettings, simulinkSettings,...
            LO,KF,AKF,DKF,GDF,...
            Cd, configuredExternally] = configure(1,"simulate",false,"Ssteel");
    plotSettings.plotBodePeak = true;
    
    % Get standard model
    MAIN;
    
    % Get identification files
    % files = {['01_100_G5.mat'],['25_1000_G100_SHORT.mat'],['250_1000_G100_SHORT'],['250_1000_G100_LONG']};
    % files = {['01_100_G5.mat'],['25_1000_G100_SHORT.mat'],['7KHz_250_1000_G100.mat'],['250_1000_G100_HighSampling'],['250_1000_G100_HighSampling']};
    files = {['1_10_G01_NEW.mat'],['100_300_G100.mat'],['100_300_G100_2.mat']};
    
    
    % Get Identified bode plot
    figure();
    
    for i = 1:length(files)
        load(files{i})
        T{i} = data{1}.Values.Time;
        dT{i} = diff(T{i});
        if std(dT{i}) > 1e-5
            warning('Sampling time not constant!')
        end
    
        sampleFreq = 1e4;
    
        if i == 1
            matchingGain = 1/db2mag(-48.6);
            minFreq = 1;
            maxFreq = 9;
            identificationGain = 0.1;
            temp = data{1}.Values.Data;
            Input = double(squeeze(temp))/matchingGain;
            Output = squeeze(data{2}.Values.Data);
            DAT = iddata(Output,Input,1/mean(dT{i}));
            nSamples = 1e3;

        elseif i == 2
            matchingGain = 1/db2mag(-46);
            minFreq = 150;
            maxFreq = 250;
            identificationGain = 50;
            temp = data{1}.Values.Data;
            Input = double(squeeze(temp))/matchingGain;
            Output = squeeze(data{2}.Values.Data);
            DAT = iddata(Output,Input,1/mean(dT{i}));
            nSamples = 5e3;
            
        elseif i == 3
            matchingGain = 1/db2mag(-46);
            minFreq = 150;
            maxFreq = 250;
            identificationGain = 100;
            temp = data{1}.Values.Data;
            Input = double(squeeze(temp))/matchingGain;
            Output = squeeze(data{2}.Values.Data);
            DAT = iddata(Output,Input,1/mean(dT{i}));
            nSamples = 5e3;
            
        end
    
        samples = logspace(log10(minFreq),log10(maxFreq),nSamples)*2*pi;
        txy = tfestimate(Input,Output,length(Output),0,samples/(2*pi),sampleFreq);
        G = frd(txy,samples,1/sampleFreq);

        bodemag(G,P)
        hold on
        
        % Plot peak in mag response
        if plotSettings.plotBodePeak == true
            if i == 1
                magnitude = squeeze(abs(G.response));
                [peakMag,I] = max(magnitude);
                peakFreq = samples(I)/(2*pi);
                plot(peakFreq,mag2db(peakMag),'o','Color',[0 0.4470 0.7410])
                text(peakFreq*1.1,mag2db(peakMag),['[',num2str(round(peakFreq,2)),',',num2str(round(mag2db(peakMag),2)),']'],'HorizontalAlignment','left')
            end
        end
   
    end
    
    modelMat(1).showBode(gca,1,1);
    title('Modelled and Identified plant')
    h = legend('1-9Hz Identification','','150-200Hz Identification','150-200Hz Identification','Model','');
    set(h,'FontSize',12);
    model = modelMat(1).dsys_sim(1,1);
end

%%
% parameters for PID control
designControl = false;
if designControl == true
    fc = simulinkSettings.controlBandwidth;
    wc = fc*2*pi;
    Kp = 1/(1*norm(freqresp(model,fc,'Hz'),2));
    wd = wc/3;
    wi = wc/5e3;
    wt = wc*3;
    
    % final PID controller
    s = tf([1,0],[0,1]);
    C = Kp*(1+wi/s)*((s/wd+1)/(s/wt+1));
    
    % Discretise and save
    Cd_sim = c2d(C,simulinkSettings.PIDSampleTime,'ZOH');
    Cd_actual = c2d(C,simulinkSettings.PIDSampleTime,'ZOH');
    % Cd = Cd*lp;
    
    Cd = ss(Cd_sim);
    
    save('./../c2000/Cd',"Cd")
    
    figure()
    margin(Cd*model,P)
end

%% Functions
% Function to plot bode for specific configuration
function evaluateParameters(Ax,parameter,range,N,baseMaterial,output,input)
    if isempty(Ax)
        figure()
        hold on
        Ax = gca;
    end

    if length(parameter) == 1 && length(range) == 2
        values = linspace(range(1),range(2),N);
        cmap = colormap('jet');
        cmap = cmap(1:floor(length(cmap)/(N-1)):length(cmap),:);
        cmap = cmap;
        opacity = 0.75;
    else
        values = range;
        cmap = [0, 0.4470, 0.7410;
                0.8500 0.3250 0.0980;
                0.9290 0.6940 0.1250;
                0.4940 0.1840 0.5560;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330];
        opacity = 1;
    end

    for i = 1:N
        if length(parameter) == 1
            value = values(i);
            color = cmap(i,:);
        else
            value = values;
            lines = length(Ax.Children);
            color = cmap(1,:);
        end

        [Beam, sBeam, modelSettings,...
        simulationSettings, plotSettings,~,~,~,~,~,~,~,~] = configure(1,parameter,value,baseMaterial);

        SYS = model(); % create model object
        [SYS,~] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 

%         SYS.showBeam([],[]);
        SYS.showBode(Ax,output,input,[color,opacity]);
    end
    
    if isempty(parameter)
        parameter = "n";
    end
    switch parameter
        case "L"
            title(strjoin([string(baseMaterial),' length variation: ',num2str(range(1)*1e3),' - ',num2str(range(end)*1e3),' mm']))
        case "b"
            title(strjoin([baseMaterial,' width variation: ',num2str(range(1)*1e3),' - ',num2str(range(end)*1e3),' mm']))
        case "h"
            title(strjoin([baseMaterial,' thickness variation: ',num2str(range(1)*1e3),' - ',num2str(range(end)*1e3),' mm']))
        case "E"
            title(strjoin([baseMaterial,' E modulus variation: ',num2str(range(1)/1e9),' - ',num2str(range(end)/1e9),' Gpa']))
        case "rho"
            title(strjoin([baseMaterial,' density variation: ',num2str(range(1)),' - ',num2str(range(end)),' Kg/m3']))
        case "stageMass"
            title(strjoin([baseMaterial,' stageMass variation: ',num2str(range(1)*1e3),' - ',num2str(range(end)*1e3),' g']))
        case "patchL"
            title(strjoin([baseMaterial,' patch length variation: ',num2str(range(1)*1e3),' - ',num2str(range(end)*1e3), 'mm']))
        otherwise
            title(['Configuration comparison'])
    end
end

% Optimisation function for patch length and position
function [peakMag] = optfnc(args)
    if any(sign(args) == -1) || args(1)*100e-3+args(2) > 100e-3
        peakMag = nan;
        peakFreq = 0;
    else
        [Beam, sBeam, modelSettings,...
            simulationSettings, plotSettings,~,~,~,~,~,~,~,~] = configure(1,["patches","patchL"],args,"Ssteel");
        
        SYS = model(); % create model object
        [SYS,~] = buildBeam(modelSettings,simulationSettings,plotSettings,Beam,sBeam,SYS); % Fill model object 
        
%         SYS.showBeam([]);

        % Find maximum of second mode peak
        fmin = 190;
        fmax = 210;
        omeg = logspace(log10(fmin*2*pi),log10(fmax*2*pi),1e3);
        sysg = frd(SYS.dsys_sim(2,1),omeg);
        [pks,lcs] = findpeaks(squeeze(abs(sysg.response)));
        peakMag = -pks(1);
        peakFreq = omeg(lcs)/(2*pi);
    end
end