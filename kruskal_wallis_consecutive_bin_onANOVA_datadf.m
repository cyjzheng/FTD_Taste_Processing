% kruskal-wallis to compare a window (ks_end - ks_start) across all tastes.
%for each trialtype, responses (testV) are averaged across trials
function [sig_p_kw_bin]=kruskal_wallis_consecutive_bin_onANOVA_datadf(datadf_PSTH, Binsize, ks_start, ks_end, Pre, alpha,anova_uni_sig_bw,uni_sig_bw,varargin)
% assigned_tastes={'S','N','Q','C'};
% Binsize=250;
% ks_start=0;
% ks_end=5000;
% Pre=-6000;
% alpha=0.05;
% pv: a pv file from a best window folder. nNeurons x nTrialTypes double
% uni_sig_bw: a file from a best window folder. n x nTrialTypes  double
%varargin{1} is meet_fr_threshold; if entered, will adjust the # total units for proportion_final
%varargin{2} is sig_kw
sig_kw=[];
p_kw=[];
results_kw = [];
neuronIDs=unique(datadf_PSTH.neuron_numeric_ID);
nNeurons=length(neuronIDs);
font="Arial";
fontsz=12;
for i= 1:nNeurons
    testV=[];
    groupV=[];
    n1_df=datadf_PSTH(datadf_PSTH.neuron_numeric_ID==i,:);
    all_scmatrix=cell2mat(n1_df.(append('scmatrix_',num2str(Binsize))));
    temp=all_scmatrix./(Binsize/1000);
    testV=temp(:,(ks_start-Pre)/Binsize+1:(ks_end-Pre)/Binsize);
    groupV =[n1_df.Taste_numeric_ID];
    for k=1:size(testV,2)
        p_kw_bin(i,k) = kruskalwallis(testV(:,k), groupV, 'off');
    end
end
sig_p_kw_bin=(p_kw_bin<alpha);
[~,id]=unique(datadf_PSTH.neuron_numeric_ID);
m=table2array(datadf_PSTH(id,'mutant')); 
a = anova_uni_sig_bw(:,1);
anova_list=false(size(m));
anova_list(a)=1;
Q331K_sig_kw_bin=sum(sig_p_kw_bin(m==1&anova_list,:),1)./nnz(m&anova_list); %taste selective/taste resp for Q331K
nonTg_sig_kw_bin=sum(sig_p_kw_bin(m==0&anova_list,:),1)./nnz(m==0&anova_list);%taste selective/taste resp for nonTg
% inc_Q331K_ind=uni_sig_bw(uni_sig_bw(:,2)==1&uni_sig_bw(:,3)==0&uni_sig_bw(:,4)==1,1);
% inc_nonTg_ind=uni_sig_bw(uni_sig_bw(:,2)==1&uni_sig_bw(:,3)==0&uni_sig_bw(:,4)==0,1);
% Q331K_sig_inc_kw_bin=sum(sig_p_kw_bin(inc_Q331K_ind,:),1)./size(inc_Q331K_ind,1);
% nonTg_sig_inc_kw_bin=sum(sig_p_kw_bin(inc_nonTg_ind,:),1)./size(inc_nonTg_ind,1);
% dec_Q331K_ind=uni_sig_bw(uni_sig_bw(:,2)==0&uni_sig_bw(:,3)==1&uni_sig_bw(:,4)==1,1);
% dec_nonTg_ind=uni_sig_bw(uni_sig_bw(:,2)==0&uni_sig_bw(:,3)==1&uni_sig_bw(:,4)==0,1);
% Q331K_sig_dec_kw_bin=sum(sig_p_kw_bin(dec_Q331K_ind,:),1)./size(dec_Q331K_ind,1);
% nonTg_sig_dec_kw_bin=sum(sig_p_kw_bin(dec_nonTg_ind,:),1)./size(dec_nonTg_ind,1);
numbins=(ks_end-ks_start)/Binsize; %to separate 0-6s evenly in several time bins
bins = linspace(ks_start,ks_end,numbins+1);
timepoint=bins(1:end-1)/1000+Binsize/1000/2;
colors=distinguishable_colors(8);

% h=gobjects(1,3);
% p=tiledlayout(1,3,"TileSpacing","compact");
% h(1,1)=nexttile;
% plot(timepoint,Q331K_sig_kw_bin,'LineWidth',2,'Color',colors(5,:));
% hold on
% plot(timepoint,nonTg_sig_kw_bin,'LineWidth',2,'Color',colors(4,:));
% title('All Taste-Selective','FontSize',fontsz,'FontName',font)
% pbaspect([1 1 1])
% h(1,2)=nexttile;
% plot(timepoint,Q331K_sig_inc_kw_bin,'LineWidth',2,'Color',colors(5,:));
% hold on
% plot(timepoint,nonTg_sig_inc_kw_bin,'LineWidth',2,'Color',colors(4,:));
% title('Taste-Selective Increase','FontSize',fontsz,'FontName',font)
% pbaspect([1 1 1])
% h(1,3)=nexttile;
% plot(timepoint,Q331K_sig_dec_kw_bin,'LineWidth',2,'Color',colors(5,:));
% hold on
% plot(timepoint,nonTg_sig_dec_kw_bin,'LineWidth',2,'Color',colors(4,:));
% title('Taste-Selective Decrease','FontSize',fontsz,'FontName',font)
% pbaspect([1 1 1])
% xlabel(p,'Time(s)')
% ylabel(p,'Porportion of responses')
% title(p,'Temporal Distribution of Selective Neurons','FontSize',fontsz,'FontName',font)
% xlim(h,[0,ks_end/1000])
% linkaxes(h,'xy');
% set(h,'TickDir','out','color','none','box','off');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% figurename1='temporal_kw_selective_on_ANOVA_inc_dec';
% saveas(gcf,figurename1);
% exportgraphics(gcf,[figurename1,'.pdf'],ContentType="vector");
% save(['temporal_sig_p_kw_on_ANOVA_bin',num2str(Binsize),'_ks_end',num2str(ks_end)],"sig_p_kw_bin");


plot(timepoint,Q331K_sig_kw_bin,'LineWidth',2,'Color',colors(5,:));
hold on
plot(timepoint,nonTg_sig_kw_bin,'LineWidth',2,'Color',colors(4,:));
title('All Taste-Selective','FontSize',fontsz,'FontName',font)
pbaspect([1 1 1])
set(gca,'TickDir','out','color','none','box','off');
set(gcf,'units','normalized','outerposition',[0 0 1 1])


%%%%%%%%%%%%%%%ANOVA%%%%%%%%%%%%%%
  data=horzcat(Q331K_sig_kw_bin,nonTg_sig_kw_bin)'; % avgFRmatrix: 1 x 4 cell (Q331K increase, Q decrease, nonTg imcrease, nonTg decrease) each with # of neurons x #of timebins double
    times=repmat(timepoint,1,2)';
    g_Q331K =repmat({'Q331K'},size(Q331K_sig_kw_bin,1)*size(Q331K_sig_kw_bin,2),1);
    g_nonTg =repmat({'nonTg'},size(nonTg_sig_kw_bin,1)*size(nonTg_sig_kw_bin,2),1);
    genotypes=char(vertcat(g_Q331K,g_nonTg));
    [pval,tbl,stats] = anovan(data,{times genotypes},'Model','full','varnames',{'times', 'genotypes'},'display','on');


end