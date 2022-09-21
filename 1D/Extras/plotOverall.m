clear
close all
set(0,'defaultTextInterpreter','latex'); 
addpath('dataSets')

load('modelMat_FINALCL1_1')
systems(1) = modelMat(1);
clear modelMat
load('modelMat_FINALCL2')
systems(2) = modelMat(1);
clear modelMat
load('modelMat_FINALOL1')
systems(3) = modelMat(1);
clear modelMat
load('modelMat_FINALOL2')
systems(4) = modelMat(1);
clear modelMat

references = cell(length(systems),1);
trues = cell(length(systems),1);
AKFest = cell(length(systems),1);
DKFest = cell(length(systems),1);
GDFest = cell(length(systems),1);

limsFull = [5,53;
          4.43,57.7;
          3.5,49;
          3.75,49];
limsLeft = [13.9,15.5;
          10.6,11.5;
          9.8,10.3;
          9.8,10.3];
limsRight = [15.01,15.04;
            10.65,10.75;
            9.99,10.02;
            9.99,10.02];

figure()
hold on
sgtitle('Overal performance')
N = 3;
for i = 1:length(systems)
    load(systems(i).simulationSettings.dataset);
    trues{i} = systems(i).simulationData.yfull(1,:);

    T = 1e4;
    startSample = (limsFull(i,1)-limsFull(i,1)+1)*T;
    stopSample = (limsFull(i,2)-limsFull(i,1)+1)*T;

    t = systems(i).simulationData.t(startSample:stopSample);

    references{i} = data{3}.Values.Data';
    AKFest{i} = systems(i).simulationData.yfull_AKF(1,:);
    DKFest{i} = systems(i).simulationData.yfull_DKF(1,:);
    GDFest{i} = systems(i).simulationData.yfull_GDF(1,:);

    subplot(length(systems),N,(1:N-1)+(i-1)*N)
        hold on
        ylabel('[m]')
        title(['Dataset ',num2str(i)])
        xlim([limsLeft(i,1),limsLeft(i,2)])
        plot(t,trues{i}(startSample:stopSample),'k')
        plot(t,AKFest{i}(startSample:stopSample),'Color',[0.4940 0.1840 0.5560])
        plot(t,DKFest{i}(startSample:stopSample),'Color',[0.3010 0.7450 0.9330]);
        plot(t,GDFest{i}(startSample:stopSample),'Color',[0.6350 0.0780 0.1840]);
        xline(limsRight(i,1),'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
        xline(limsRight(i,2),'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)

        if i == length(systems)
            xlabel('time [s]')
        end

        if i == 1
            legend('True','AKF','DKF','GDF','','','Location','northwest')
        end

    subplot(length(systems),N,i*N)
        hold on
        ylabel('[m]')
        xlim([limsRight(i,1),limsRight(i,2)])
        plot(t,trues{i}(startSample:stopSample),'k')
        plot(t,AKFest{i}(startSample:stopSample),'Color',[0.4940 0.1840 0.5560])
        plot(t,DKFest{i}(startSample:stopSample),'Color',[0.3010 0.7450 0.9330]);
        plot(t,GDFest{i}(startSample:stopSample),'Color',[0.6350 0.0780 0.1840]);
        if i == length(systems)
            xlabel('time [s]')
        end

end

means = zeros(2,1);
meansP = zeros(2,1);

%%
figure() % Input estimate figure!
hold on
sgtitle('Input estimation performance')
N = 3;

LimsFull = [3.5,49];
LimsLeft = [9,11];
LimsRight = [9.9975,10.01];
matchingGain = 1/db2mag(-48.6);

startSample = (limsFull(i,1)-limsFull(i,1)+1)*T;
stopSample = (limsFull(i,2)-limsFull(i,1)+1)*T;

for i = 3:4
    load(systems(i).simulationSettings.dataset);
    AKFinEst = systems(i).simulationData.qfull_AKF(end,:);
    DKFinEst = systems(i).simulationData.ufull_DKF;
    GDFinEst = systems(i).simulationData.ufull_GDF;
    trueIn = data{4}.Values.Data'/(matchingGain*1e3);
    t = systems(i).simulationData.t(startSample:stopSample);
    
    j = i-2;
    subplot(2,N,(1:N-1)+(j-1)*N)
        hold on
        ylabel('Estimated force [N]')
        title(['Dataset ',num2str(i)])
        xlim([limsLeft(i,1),limsLeft(i,2)])
        plot(t,AKFinEst(startSample:stopSample),'Color',[0.4940 0.1840 0.5560])
        plot(t,DKFinEst(startSample:stopSample),'Color',[0.3010 0.7450 0.9330]);
        plot(t,GDFinEst(startSample:stopSample),'Color',[0.6350 0.0780 0.1840]);
        xline(limsRight(i,1),'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
        xline(limsRight(i,2),'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
        if j == 1
            legend('AKF','DKF','GDF','','','Location','northwest')
        end
        if i == length(systems)
            xlabel('time [s]')
        end
    subplot(2,N,j*N)
        hold on
        grid on
        ylabel('[m]')
        xlim([limsRight(i,1),limsRight(i,2)])
        plot(t,AKFinEst(startSample:stopSample),'Color',[0.4940 0.1840 0.5560])
        plot(t,DKFinEst(startSample:stopSample),'Color',[0.3010 0.7450 0.9330]);
        plot(t,GDFinEst(startSample:stopSample),'Color',[0.6350 0.0780 0.1840]);
        if i == length(systems)
            xlabel('time [s]')
        end
end
