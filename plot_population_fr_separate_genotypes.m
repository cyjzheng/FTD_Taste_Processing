function plot_population_fr_separate_genotypes(xlim_start,xlim_finish,Pre,Post,anova_start,anova_end,Binsize,FRMatrix,assigned_tastes)
fontSize=12;
font='Arial';
sm_fr_all={};
ste_fr_all={};

for plotn = 1 : 4
    frr=[];
    reshaped_FRMatrix =cat(3, FRMatrix{plotn,:});
frr=mean(reshaped_FRMatrix,3);
   averaged_FRMatrix{plotn}= frr;
     for k=1:size( frr,1)
             frr(k,:)=gaussmooth( frr(k,:),1000,2);
     end
     fr=mean(frr,1);
     ste_fr=std(frr,0,1)/sqrt(size(frr,1));
     sm_fr_all{1,plotn}=fr;
     ste_fr_all{1,plotn}=ste_fr;
end
% plot
Q331K_fr_inc=vertcat(sm_fr_all{:,1});
Q331K_ste_inc=vertcat(ste_fr_all{:,1});
nonTg_fr_inc=vertcat(sm_fr_all{:,3});
nonTg_ste_inc=vertcat(ste_fr_all{:,3});
std_option=1;
plot_firing_rate_two_genotypes (Pre, Post, Binsize, Q331K_fr_inc,Q331K_ste_inc, nonTg_fr_inc, nonTg_ste_inc, fontSize,font, xlim_start,xlim_finish, assigned_tastes,std_option )
f=gcf;
f.Position=get(0,'ScreenSize');
figurename1=['normalized_avg_separate_genotype_increase_',num2str(Binsize),'ms'];
savefig(gcf,figurename1);
exportgraphics(gcf,[figurename1,'.pdf'],ContentType="vector");
Q331K_fr_dec=vertcat(sm_fr_all{:,2});
Q331K_ste_dec=vertcat(ste_fr_all{:,2});
nonTg_fr_dec=vertcat(sm_fr_all{:,4});
nonTg_ste_dec=vertcat(ste_fr_all{:,4});
plot_firing_rate_two_genotypes (Pre, Post, Binsize, Q331K_fr_dec,Q331K_ste_dec, nonTg_fr_dec, nonTg_ste_dec, fontSize,font, xlim_start,xlim_finish, assigned_tastes,std_option )
figurename2=['normalized_avg_separate_genotype_decrease_',num2str(Binsize),'ms'];
f=gcf;
f.Position=get(0,'ScreenSize');
savefig(gcf,figurename2);
exportgraphics(gcf,[figurename2,'.pdf'],ContentType="vector");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%inc/dec ratio%%%%%%%%%%%%%%%%%%%%%%%%%
colors=distinguishable_colors(8);
Q331K_EI_diff=abs(Q331K_fr_inc-0.5)-abs(Q331K_fr_dec-0.5);
nonTg_EI_diff=abs(nonTg_fr_inc-0.5)-abs(nonTg_fr_dec-0.5);
timept=linspace(Pre/1000,Post/1000,(Post-Pre)/Binsize+1);
timepoint=timept(1:end-1)+Binsize/1000/2;
figure;
plot(timepoint,Q331K_EI_diff,'w','LineWidth',1,Color=colors(5,:));hold on;
plot(timepoint,nonTg_EI_diff,'w','LineWidth',1,Color=colors(4,:))
    xlim([xlim_start/1000 (xlim_finish)/1000]);
 ylabel('E/I Difference','FontSize', fontSize,'FontName',font);

    xlabel('Time(s)','FontSize', fontSize,'FontName',font)
  pbaspect([1 1 1]);
    set(gca,'TickDir','out','color','none','box','off');
figurename3=['EI_difference',num2str(Binsize),'ms'];
savefig(gcf,figurename3);
exportgraphics(gcf,[figurename3,'.pdf'],ContentType="vector");
% Q331K_EI_ratio=gaussmooth(abs(Q331K_fr_inc-0.5)./abs(0.5-Q331K_fr_dec),1000,2);
% nonTg_EI_ratio=gaussmooth(abs(nonTg_fr_inc-0.5)./abs(0.5-nonTg_fr_dec),1000,2);
% plot(timepoint,Q331K_EI_ratio,'w','LineWidth',1,Color=colors(5,:));hold on;
% plot(timepoint,nonTg_EI_ratio,'w','LineWidth',1,Color=colors(4,:))
%   xlim([xlim_start/1000 (xlim_finish)/1000]);
%  ylabel('E/I Ratio','FontSize', fontSize,'FontName',font);
% 
%     xlabel('Time(s)','FontSize', fontSize,'FontName',font)
%   pbaspect([1 1 1]);
%     set(gca,'TickDir','out','color','none','box','off');
% figurename3=['EI_ratio_smoothed',num2str(Binsize),'ms'];
% savefig(gcf,figurename3);
% exportgraphics(gcf,[figurename3,'.pdf'],ContentType="vector");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ANOVA for each genotype%%%%%%%%%%%%%%%%%%%%%%%%%
trialtype=assigned_tastes;
Timebin_Labels=num2cell((anova_start-Pre)/Binsize+1:(anova_end-Pre)/Binsize);
%for excitatory responses
for i=1:size(FRMatrix,1)
        anovaFRmatrix{i}=averaged_FRMatrix{i}(:,(anova_start-Pre)/Binsize+1:(anova_end-Pre)/Binsize);
        nNeuron=size(FRMatrix{i},1);
        Timebin_temp=repmat(Timebin_Labels,nNeuron,1);
        Timebin_m{i}=reshape(Timebin_temp,[],1);
end

inc_Q331K=reshape(vertcat(anovaFRmatrix{1}),[],1);inc_nonTg=reshape(vertcat(anovaFRmatrix{3}),[],1);
data=vertcat(inc_Q331K,inc_nonTg);
    times=cell2mat(vertcat(vertcat(Timebin_m{1}),vertcat(Timebin_m{3})));
    g_Q331K =repmat({'Q331K'},size(inc_Q331K,1)*size(inc_Q331K,2),1);
    g_nonTg =repmat({'nonTg'},size(inc_nonTg,1)*size(inc_nonTg,2),1);
    genotypes=char(vertcat(g_Q331K,g_nonTg));
    [pval,tbl,stats] = anovan(data,{times genotypes},'Model','full','varnames',{'times', 'genotypes'},'display','on');
% Convert table to string
    tblStr = formattedDisplayText(tbl);
   fid = fopen(append('anovaTable_inc',num2str(anova_start), '_', num2str(anova_end),'.txt'), 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
    save(append('ANOVA_inc_',num2str(anova_start), '_', num2str(anova_end)),'stats','pval');


dec_Q331K=reshape(vertcat(anovaFRmatrix{2}),[],1);dec_nonTg=reshape(vertcat(anovaFRmatrix{4}),[],1);
data=vertcat(dec_Q331K,dec_nonTg);
    times=cell2mat(vertcat(vertcat(Timebin_m{2}),vertcat(Timebin_m{4})));
    g_Q331K =repmat({'Q331K'},size(dec_Q331K,1)*size(dec_Q331K,2),1);
    g_nonTg =repmat({'nonTg'},size(dec_nonTg,1)*size(dec_nonTg,2),1);
    genotypes=char(vertcat(g_Q331K,g_nonTg));
    [pval,tbl,stats] = anovan(data,{times genotypes},'Model','full','varnames',{'times', 'genotypes'},'display','on');
% Convert table to string
    tblStr = formattedDisplayText(tbl);
   fid = fopen(append('anovaTable_dec',num2str(anova_start), '_', num2str(anova_end),'.txt'), 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
    save(append('ANOVA_dec_',num2str(anova_start), '_', num2str(anova_end)),'stats','pval');

 
end