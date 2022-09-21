clear all;
close all;
addpath("dataSets")

load('modelMat_FINALCL1_1')
sys1 = modelMat(1);
clear modelMat
load('modelMat_FINALCL2')
sys2 = modelMat(1);
clear modelMat

[meanAccelerationSections1,meanErrSections1] = sys1.plotFinalResults();
[meanAccelerationSections2,meanErrSections2] = sys2.plotFinalResults();

mAc = [meanAccelerationSections1, meanAccelerationSections2];
mEr = [meanErrSections1, meanErrSections2];


figure()
    hold on
    sgtitle('Mean absolute acceleration and filter error values')
    subplot(3,1,1)
        hold on
        bar(mAc')
        xlim([0.5,6.5])
    subplot(3,1,[2 3])    
        hold on
        b = bar(mEr','FaceColor','flat');
        b.CData(2,:) = [0.4940 0.1840 0.5560];
        xlim([0.5,6.5])
        xticklabels([1:6])
        ylabel('')
        xlabel('Data set section')




AKFcorrCoeff = corrcoef(mEr(1,:)',mAc');
DKFcorrCoeff = corrcoef(mEr(2,:)',mAc');
GDFcorrCoeff = corrcoef(mEr(3,:)',mAc');

figure();
    hold on
    title('Acceleration versus filter error')
    xlabel('Mean absolute acceleration signal value [m/s$^2$]')
    ylabel('Mean absolute filter error [m]')
    
    plot(mAc,mEr(1,:),'o','Color',[0.4940 0.1840 0.5560]);
    plot(mAc,mEr(2,:),'d','Color',[0.3010 0.7450 0.9330]);
    plot(mAc,mEr(3,:),'^','Color',[0.6350 0.0780 0.1840]);
    
    xlimits = xlim();
    pAKF = polyfit(mAc,mEr(1,:),1);
    fAKF = polyval(pAKF,mAc);
    
    pDKF = polyfit(mAc,mEr(2,:),1);
    fDKF = polyval(pDKF,mAc);
    
    pGDF = polyfit(mAc,mEr(3,:),1);
    fGDF = polyval(pGDF,mAc);
    
    plot(mAc,fAKF,'Color',[0.4940 0.1840 0.5560]);
    plot(mAc,fDKF,'Color',[0.3010 0.7450 0.9330]);
    plot(mAc,fGDF,'Color',[0.6350 0.0780 0.1840]);
    
    legend('AKF','DKF','GDF','','Location','northwest');


