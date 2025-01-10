load('data_df_PSTH');
assigned_tastes={'S','N','Q','C'};
stimulus='valve_on';
window=250;
pre_window=1000;
post_window=250;
Post=10000;
Pre=-6000;
step=25;
limit=5000;
alpha=0.05;
n_sd=2; % set criterion for how good the best window is. response is n_sd*standard deviation compared to the baseline
foldername='202210_4taste_pre6000_post10000';
datasavepath='D:\AF Lab\matlab_data\best_window';
subfoldername=[foldername,'_final'];
best_window_final(assigned_tastes,stimulus,window,pre_window,post_window,Pre,step, limit,alpha,n_sd,datadf_PSTH,datasavepath,subfoldername);
Binsize=250;
baseline_window=2000; %baseline window
xlim_finish=8000;
xlim_start=-4000;
baseline_option=1;
anova_start=0;
anova_end=5000;
frsavepath='D:\AF Lab\matlab_data\fr_plots';
cd([datasavepath,'\',subfoldername,'\SD',num2str(n_sd),'\bw_avgALL_responsive_',stimulus,'_pre',num2str(pre_window),'_post',num2str(window),'_postwindow',num2str(post_window),'_limit',num2str(limit)])
%averaged increase/decrease firing
%plot mean proportion across all tastes
load('proportional_final');
res=reshape(mean(proportion_final,1),[],2)';
groupLabels={'Q331K','non-Tg'};
b= bar(res,'hist');
set(gca,'XTickLabel',groupLabels);
colors=linspecer(9);
ccl{1,1}=colors(9,:);
ccl{2,1}=colors(8,:);
set(b,{'FaceColor'},ccl)
ylabel('Proportion Taste Responsive');
b(1,1).DisplayName='Inc. Firing';
b(1,2).DisplayName='Dec. Firing';
legend(b,'location','eastoutside');
axis square
set(gca,'TickDir','out','color','none','box','off');
figurename=['average_response_final_proportion.pdf'];
exportgraphics(gca,figurename,'BackgroundColor','none');
load('uni_sig_bw');
norm=1;
plot_population_avg_fr_rate_bw_cz(uni_sig_bw,n_sd,assigned_tastes,datadf_PSTH,stimulus,frsavepath,foldername,baseline_window,Binsize,Pre,Post,norm,xlim_finish,xlim_start,baseline_option,anova_start, anova_end,pre_window);
load('FRMatrix_per_plot_taste');
plot_population_fr_separate_genotypes(xlim_start,xlim_finish,Pre,Post,anova_start,anova_end,Binsize,FRMatrix,assigned_tastes)
