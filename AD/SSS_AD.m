
load('AD_SSS.mat')
SS_AD = SSS;
load('CN_SSS.mat');
SS_CN = SSS;
load('LMCI_SSS.mat');
SS_LMCI = SSS;
SSS = [SS_AD;SS_LMCI;SS_CN];
x=1:9;
bar(x,SSS)
ylim([0.0001,1000])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'XTick',1:9,'xticklabel',{'\lambda_{A\beta}','\lambda_\tau'...
    ,'\lambda_{N\tau_p}','\lambda_{CN}','\lambda_{C\tau}',...
    'K_{A\beta}','K_{\tau_p}','K_N','K_C'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
lgd=legend('AD','LMCI','CN');
lgd.Box='off';
lgd.FontSize = 15;
lgd.ItemTokenSize = [10,6];
box off
