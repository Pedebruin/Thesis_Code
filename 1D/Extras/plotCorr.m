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

for i = [1,4]
figure()
    title('Filter error and Mean absolute acceleration')
    hold on

    yyaxis left
    br = bar([i:i+2],mAc(i:i+2)');
    br.FaceAlpha = 0.2;
    br.EdgeAlpha = br.FaceAlpha;
    br.EdgeColor = [0 0.4470 0.7410];
%     br.ZData = ones(size(br.XData));


    ylim([0,max(mAc)])
    ylabel('mean absolute acceleration [m/s$^2$]')

    yyaxis right
    bl = bar([i:i+2],mEr(:,i:i+2)','FaceColor','flat');
    for j = 1:3
        bl(1).CData(j,:) = [0.4940 0.1840 0.5560];
        bl(2).CData(j,:) = [0.3010 0.7450 0.9330];
        bl(3).CData(j,:) = [0.6350 0.0780 0.1840];

        bl(1).EdgeColor = [0.8500 0.3250 0.0980];
        bl(2).EdgeColor = [0.8500 0.3250 0.0980];
        bl(3).EdgeColor = [0.8500 0.3250 0.0980];

%         bl(1).ZData = zeros(size(bl(1).XData));
%         bl(2).ZData = zeros(size(bl(2).XData));
%         bl(3).ZData = zeros(size(bl(3).XData));
    end

    xlim([i-0.5,i+3-0.5])
    ylim([0,max(max(mEr))*1.5])

    xticks(i:i+2)
    ylabel('')
    xlabel('Data set section')
    ylabel('mean absolute filter error [m]')
    grid on

    set(gca, 'SortMethod', 'depth')
end

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


