load('PSTH_auROC_bin250.mat');
load('mutant_info.mat'); %nNeurons+ x 1 double. indicate mutant =1; nonTg=0;
Binsize=250;
assigned_tastes={'S','N','Q','C'};
Pre=-6000;
Post=10000;
limit=5000;
window=250;
pre_window=4000; %how much time before valve on is plotted
baseline_window=4000;
nNeurons=size(PSTH_auROC,1);
PI_ind=zeros(nNeurons,1);
PI_positive_ind=zeros(nNeurons,1);
 PI_positive_ind_bin=zeros(nNeurons, size(-Pre/Binsize+1:(-Pre+limit)/Binsize,2));
for i= 1:nNeurons
    all_TrialType_auROC=squeeze(PSTH_auROC(i,:,:));
    LR_same=0.5*(abs(log(all_TrialType_auROC(1,:)./all_TrialType_auROC(2,:)))+abs(log(all_TrialType_auROC(3,:)./all_TrialType_auROC(4,:))));
    LR_oppo=0.25*(abs(log(all_TrialType_auROC(1,:)./all_TrialType_auROC(4,:)))+abs(log(all_TrialType_auROC(1,:)./all_TrialType_auROC(3,:)))+abs(log(all_TrialType_auROC(2,:)./all_TrialType_auROC(4,:)))+abs(log(all_TrialType_auROC(2,:)./all_TrialType_auROC(3,:))));
    PI_all(i,:)=LR_oppo-LR_same;
    baseline_avg=mean(PI_all(i,(-Pre-4000)/Binsize+1:(-Pre-4000+baseline_window)/Binsize));
    baseline_std=std(PI_all(i,(-Pre-4000)/Binsize+1:(-Pre-4000+baseline_window)/Binsize),[],"all");
    PI_ind(i,1)=any(abs(PI_all(i,-Pre/Binsize+1:(-Pre+limit)/Binsize))>6*baseline_std);
    PI_positive_ind(i,1)=any((PI_all(i,-Pre/Binsize+1:(-Pre+limit)/Binsize)-abs(baseline_avg))>6*baseline_std);
 PI_positive_ind_bin(i,:)=(PI_all(i,-Pre/Binsize+1:(-Pre+limit)/Binsize)-abs(baseline_avg))>6*baseline_std; % indicates when the neuron has significant and positive PI, during 0 -5 s
end
save("PI_all_ind","PI_all","PI_ind","PI_positive_ind_bin","PI_positive_ind");
%add taste responsive information
cd('D:\AF Lab\matlab_data\best_window\202210_4taste_pre6000_post10000_goodcentral_vo53094_v2b1a\SD2\bw_avgALL_responsive_valve_on_pre1000_post250_postwindow250_limit5000')
load('uni_sig_bw.mat')
m(uni_sig_bw(:,1),2)=1;
cd('D:\AF Lab\matlab_data\best_window\202210_4taste_pre6000_post10000_goodcentral_vo53094_v2b1a\SD2\bw_avgALL_responsive_valve_on_pre1000_post250_postwindow250_limit5000\double_bw_responsive_alltrialtypes_first_lick_baseline1000_post500_bin25')
load('uni_sig_only_taste.mat')
m(uni_sig_only_taste(:,1),3)=1;
cd('D:\AF Lab\matlab_data\best_window\202210_4taste_pre6000_post10000_goodcentral_vo53094_v2b1a\SD2\bw_avgALL_responsive_valve_on_pre1000_post250_postwindow250_limit5000\kw_bin1000_start0_end3000')
load('kw_uni_sig_bw.mat')
%create a second column in m that indicates the unit's taste responsiveness
m(kw_uni_sig_bw(:,1),4)=1;%create a second column in m that indicates the unit's taste selectiveness 0-3s
clear kw_uni_sig_bw
cd('D:\AF Lab\matlab_data\best_window\202210_4taste_pre6000_post10000_goodcentral_vo53094_v2b1a\SD2\bw_avgALL_responsive_valve_on_pre1000_post250_postwindow250_limit5000\kw_bin1000_start3000_end5000')
load('kw_uni_sig_bw.mat');
m(kw_uni_sig_bw(:,1),5)=1;%create a second column in m that indicates the unit's taste selectiveness 3-5s
cd('D:\AF Lab\matlab_data\best_window\202210_4taste_pre6000_post10000_goodcentral_vo53094_v2b1a\SD2\bw_avgALL_responsive_valve_on_pre1000_post250_postwindow250_limit5000\anova_bin1000');
load('anova_uni_sig_bw.mat')
m(anova_uni_sig_bw(:,1),6)=1;
crit={};
crit_sig={};
sorted_crit={};
sorted_crit_sig={};
crit{1,1}=PI_all(m(:,1)==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
crit{1,2}=PI_all(m(:,1)==0,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
crit_sig{1,1}=PI_all(m(:,1)==1&PI_positive_ind==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
crit_sig{1,2}=PI_all(m(:,1)==0&PI_positive_ind==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
for k=2:size(m,2)
    crit{k,1}=PI_all(m(:,1)==1&m(:,k)==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
    crit{k,2}=PI_all(m(:,1)==0&m(:,k)==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
    crit_sig{k,1}=PI_all(m(:,1)==1&m(:,k)==1&PI_positive_ind==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
    crit_sig{k,2}=PI_all(m(:,1)==0&m(:,k)==1&PI_positive_ind==1,(-Pre-pre_window)/Binsize+1:(-Pre+limit)/Binsize);
end
for k=1:size(m,2)
    [~, index1] = sort(max(crit{k,1}, [], 2),'descend');
    sorted_crit{k,1} = crit{k,1}(index1, :);
    [maxv, index2] = sort(max(crit_sig{k,1}, [], 2),'descend');
    sorted_crit_sig{k,1} = crit_sig{k,1}(index2, :);
    [~, index3] = sort(max(crit{k,2}, [], 2),'descend');
    sorted_crit{k,2} = crit{k,2}(index3, :);
    [~, index4] = sort(max(crit_sig{k,2}, [], 2),'descend');
    sorted_crit_sig{k,2} = crit_sig{k,2}(index4, :);
    max_c(k)=max(vertcat(crit{k,1},crit{k,2}),[],'all');
    min_c(k)=min(vertcat(crit{k,1},crit{k,2}),[],'all');
    max_s(k)=max(vertcat(crit_sig{k,1},crit_sig{k,2}),[],'all');
    min_s(k)=min(vertcat(crit_sig{k,1},crit_sig{k,2}),[],'all');
end

Xlim=Pre/1000+Binsize/2/1000:Binsize/1000:Post/1000;
Xlim_of_interest=(-pre_window)/1000+Binsize/2/1000:Binsize/1000:limit/1000;

p=tiledlayout(6,2,"TileSpacing","compact");
h=gobjects(6,2);
% heatmap_variable=sorted_crit;
heatmap_variable=sorted_crit_sig;
% for k=4:size(m,2)
for k=6
    h(k,1)=nexttile;
    imagesc(Xlim_of_interest, 1, heatmap_variable{k,1});
    ylim2use = get(gca, 'ylim');
    line([0 0], ylim2use, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1)
    h(k,2)=nexttile;
    imagesc(Xlim_of_interest, 1, heatmap_variable{k,2});
    % set(h(k,:), 'Colormap', turbo, 'CLim', [min_s(k) max_s(k)])
    set(h(k,:), 'Colormap', turbo, 'CLim', [-0.1 1])
    cbh = colorbar(h(k,2));
    ylim2use = get(gca, 'ylim');
    line([0 0], ylim2use, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1)
end
exportgraphics(gcf,['Q331KnonTg_heatmap_ANOVA.pdf'])


%%%%plot proportion over taste selective units
fontSize=12;
font='Arial';
numbins=(limit-0)/Binsize; %to separate 0-6s evenly in several time bins
bins = linspace(0,limit,numbins+1);
timepoint=bins(1:end-1)/1000+Binsize/1000/2;
Q331K_anova_ind=anova_uni_sig_bw(anova_uni_sig_bw(:,4)==1,1);
nonTg_anova_ind=anova_uni_sig_bw(anova_uni_sig_bw(:,4)==0,1);
Q331K_sig_inc_kw_bin=sum(PI_positive_ind_bin(Q331K_anova_ind,:),1)./size(Q331K_anova_ind,1);
nonTg_sig_inc_kw_bin=sum(PI_positive_ind_bin(nonTg_anova_ind,:),1)./size(nonTg_anova_ind,1);
colors=distinguishable_colors(8);
b1=bar(timepoint,nonTg_sig_inc_kw_bin,'EdgeColor',colors(4,:),'EdgeAlpha',0.5,'FaceColor',colors(4,:),'FaceAlpha',0.5,'BarWidth',1);
b2=bar(timepoint,Q331K_sig_inc_kw_bin,'EdgeColor',colors(5,:),'EdgeAlpha',0.5,'FaceColor',colors(5,:),'FaceAlpha',0.5,'BarWidth',1); hold on;
xlim([0 limit/1000]);
legend([b2, b1],{'Q331K','non-Tg'},"Box","off",Location="northeast")
axis square
ylabel('Proportion of Palatability-Coding Units')
xlabel('Time (s)')
set(gca,'FontSize', fontSize,'FontName',font);
set(gca,'TickDir','out','color','none','box','off');
figurename=['Proportion_palatability_coding_',num2str(limit/1000),'s'];
f=gcf;
f.Position=get(0,'ScreenSize');
savefig(gcf,figurename);
exportgraphics(gcf,[figurename,'.pdf'],ContentType="vector");
%ranksum for time points 0-5 s after taste delivery
[p,st]=ranksum(Q331K_sig_inc_kw_bin,nonTg_sig_inc_kw_bin) ;
 
% for k=1:size(m,2)
%     Q331K_P(k,:)=mean(crit{k,1},1);
%     nonTg_P(k,:)=mean(crit{k,2},1);
%     sterr1(k,:)=std(crit{k,1},[],1)./sqrt(size(crit{k,1},1));
%     sterr2(k,:)=std(crit{k,2},[],1)./sqrt(size(crit{k,2},1));
%     yy(k,:) = [Q331K_P(k,:) + sterr1(k,:),fliplr(Q331K_P(k,:) - sterr1(k,:))];
%     yy2(k,:) = [nonTg_P(k,:) + sterr1(k,:),fliplr(nonTg_P(k,:) - sterr1(k,:))];
% end
for k=1:size(m,2)
    Q331K_P(k,:)=mean(crit_sig{k,1},1);
    nonTg_P(k,:)=mean(crit_sig{k,2},1);
    sterr1(k,:)=std(crit_sig{k,1},[],1)./sqrt(size(crit_sig{k,1},1));
    sterr2(k,:)=std(crit_sig{k,2},[],1)./sqrt(size(crit_sig{k,2},1));
    yy(k,:) = [Q331K_P(k,:) + sterr1(k,:),fliplr(Q331K_P(k,:) - sterr1(k,:))];
    yy2(k,:) = [nonTg_P(k,:) + sterr1(k,:),fliplr(nonTg_P(k,:) - sterr1(k,:))];
end
%making a shadow with standard deviation
xx = [Xlim_of_interest(1:end), fliplr(Xlim_of_interest(1:end))];

colors=distinguishable_colors(8);


fontSize = 12;
t=tiledlayout(2,3);
subPlotNames={'All units','Taste responsive', 'Taste responsive- taste only' , 'kw 0-3 s' ,'kw 3-5 s', 'ANOVA'};
% subPlotNames={['increased firing Q331K (n=',num2str(size(list{1,1},1)),')'],['decreased firing Q331K (n=',num2str(size(list{1,2},1)),')'] ,['increased firing non-Tg (n=',num2str(size(list{1,3},1)),')'] ,['decreased firing non-Tg (n=',num2str(size(list{1,4},1)),')']};
for plotn=1:size(m,2)
    ax(plotn)=nexttile;
    p(plotn,1)= fill(xx, yy(plotn,:), colors(5,:),'LineStyle','none','FaceAlpha',1); hold on;
    set(p(plotn,1),{'DisplayName'},{'Q331K'});
    p(plotn,2)=plot(Xlim_of_interest,Q331K_P(plotn,:),'w','LineWidth',1); hold on;
    set(p(plotn,2),{'DisplayName'},{''});
    p(plotn,3)= fill(xx, yy2(plotn,:), colors(4,:),'LineStyle','none','FaceAlpha',1);hold on;
    set(p(plotn,3),{'DisplayName'},{'non-Tg'});
    p(plotn,4)= plot(Xlim_of_interest,nonTg_P(plotn,:),'w','LineWidth',1); hold on;
    set(p(plotn,4),{'DisplayName'},{''});
    y=ylim;
    p(plotn,5)= plot([0 0],[y(1) y(2)],'--','Color','k');
    set(p(plotn,5),{'DisplayName'},{'valve on'});
    pbaspect([1 1 1])
    box('off')
    title(subPlotNames{plotn});
end
linkaxes(ax(:),'xy');
xlim(ax(:),[-2 limit/1000])
set(ax,'TickDir','out','color','none','box','off');
legend(ax(3),'location','northeastoutside','FontSize', fontSize);
t.TileSpacing = 'loose';
xlabel(t,'Time (s)','FontSize',12)
ylabel(t,'Δ auROC Normalized Firing Rate','FontSize',12)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
exportgraphics(gcf,['pre_window',num2str(pre_window),'_baseline_window',num2str(baseline_window),'_limit',num2str(limit),'_window',num2str(window),'_logratio_palatability_genotype_sig.pdf'])
save("crit_sig","crit_sig");
save("crit","crit");
%ANOVA for time points 0-5 s after taste delivery
Group1.data=crit_sig{6,1}(:,17:end)';
Group2.data=crit_sig{6,2}(:,17:end)'; %rows: different time points; column: different neurons within a group
Group1.name='Q331K'; % group 1 name
Group2.name='nonTg'; % group 2 name
[pval,tbl,stats] = fun_anova(Group1,Group2);
   % Convert table to string
    tblStr = formattedDisplayText(tbl);
    % Write string to file
    fid = fopen(['anovaTable_palatability.txt'], 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
    save(['ANOVA_palatability'],'stats','pval');