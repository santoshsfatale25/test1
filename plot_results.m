%% Commands to plot results for RANDOM-RANDOM policy (Simulation+Theorotical)

%% Plot hit rate Vs Producer Number for Zipf (cmp_beta)
%% Individual hit rate for cmp_beta
%% make sure that you are taking data accordingly.
nn=1:Producers;
% cd('D:\IoT\IoT\31Jan\LRU');
% load('LRU_LRU_cmp_beta'); % for modified LRU
% load('LRU_LRU_cmp_beta_plain_LRU') % for plain LRU
% load('LRU_LRU_cmp_beta_plain_LRU1') % for simple plain LRU 
for bb=1:length(FCommon);
    h=figure;
    plot(nn,hit_rate_Simul_Zipf(:,bb),'-^r','LineWidth',2);hold on;
    plot(nn,hit_rate_theor_Zipf_modi(:,bb),'-ob','LineWidth',2);hold on;
    plot(nn,hit_rate_theor_Zipf_plain(:,bb),'-*g','LineWidth',2);
    plot(nn,hit_rate_theor_Zipf_plain1(:,bb),'-*k','LineWidth',2);
    xlabel('Producer Index','Fontsize',14);
    ylabel('Hit rate','Fontsize',14);
    title(sprintf('Hit rate Vs Producer Index for FCommon=%d',FCommon(bb)),'Fontsize',12);
    legend('Hit rate (Simulation)','Hit rate modified(Theoretical)','Hit rate plain(Theoretical)','Hit rate modified1(Theoretical)');
    grid on;
    hold off;
%     saveas(h,sprintf('Hit rate Vs Producer Index for beta=%1.2f.jpg',beta(bb)));
%     saveas(h,sprintf('Hit rate Vs Producer Index for beta=%1.2f.fig',beta(bb)));
end

%%