% this version plots normalized to baseline auROC for each trial type separately
% separate responses based on inc dec and sig uni_sig (n x 4)
% removed double counted neurons (has increase and decrese responses to some trial types)
% default to smooth the curves
%inputs include: datadf_PSTH
function [FRMatrix,sm_fr_all]=plot_population_avg_fr_rate_bw_cz(uni_sig,bw_sd,assigned_tastes,datadf_PSTH,stimulus,frsavepath,foldername,baseline_window,Binsize,Pre,Post,norm,xlim_finish,xlim_start,baseline_option,anova_start, anova_end,pre_window)
% bw_sd=2; %if uni_sig is a best window result, enter the standard deviation without the decimal point. if it's not best window results, enter 0
% Binsize=250;
% Post=10000;
% Pre=-6000;
% baseline_window=2000; %baseline window
% pre_window=2000; %the criteria for processing bw
% stimulus='valve_on';
% foldername='202210_4taste_pre6000post10000goodcentral_final';
% norm=1;
% xlim_finish=8000;
% xlim_start=-4000;
% baseline_option=1;
% anova_start=0;
% anova_end=5000;
% assigned_tastes={'S','N','Q','C'};
%%
if norm==1
    figurename=['normalized_population_avg_all_',num2str(Binsize),'ms_','double_removed'];
else
    figurename=['population_avg_all_',num2str(Binsize),'ms_','double_removed'];
end
list={};
%increased firing, mutant
list{1,1}=uni_sig(uni_sig(:,2)==1&uni_sig(:,3)==0& uni_sig(:,4)==1);
%decreased firing, mutant
list{1,2}=uni_sig(uni_sig(:,2)==0&uni_sig(:,3)==1& uni_sig(:,4)==1);
%increased firing, wt
list{1,3}=uni_sig(uni_sig(:,2)==1&uni_sig(:,3)==0& uni_sig(:,4)==0);
%decreased firing, wt
list{1,4}=uni_sig(uni_sig(:,2)==0&uni_sig(:,3)==1& uni_sig(:,4)==0);
%% auROC PSTH
sm_fr_all={};
for plotn = 1 : 4
    for j=1:length(assigned_tastes)
        fr_per_taste=[];
        norm_fr_per_taste=[];
        sm_fr=[];
        fr=[];
        a=list{1,plotn};
        if ~(isempty(a))
            for n=1:length(a)
                %                 temp=getfield(sub_sum(n), (append('T_',num2str(j))),PSTHname,'FRmatrix');
                n1_df=datadf_PSTH(datadf_PSTH.neuron_numeric_ID==a(n)&strcmp(string(datadf_PSTH.assigned_tastes),assigned_tastes{j}),:);
                all_scmatrix=[cell2mat(n1_df.(append('scmatrix_',num2str(Binsize))))];
                temp=all_scmatrix./(Binsize/1000);
                if baseline_option==1 %baseline from the begining of the entire trial
                    baseline=temp(:,(-Pre-4000)/Binsize+1:(-Pre-4000+baseline_window)/Binsize); % baseline starts from -4s before taste delivery
                elseif baseline_option==2 %baseline right before stimulus
                    baseline=temp(:,(-Pre-baseline_window)/Binsize+1:(-Pre)/Binsize);
                end
                if norm==1
                    [PSTH_auROC] = psth_auROC(temp,baseline_window,Binsize,baseline);
                    norm_fr_per_taste(n,:)=PSTH_auROC;
                    FRMatrix{plotn,j}= norm_fr_per_taste;
                elseif norm==0
                    fr_per_taste(n,:)=mean(temp,1);
                    FRMatrix{plotn,j}= fr_per_taste;
                end
            end
        else
            FRMatrix{plotn,j}=NaN; % in case one category is empty
        end
        if norm==1
            fr = mean(norm_fr_per_taste,1);
        elseif norm==0
            fr = mean(fr_per_taste,1);
        end
        sm_fr= gaussmooth(fr,100,2); % smoothing the firing rate
        sm_fr_all{j,plotn}=sm_fr;
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% if you want to plot with the FRMatrix
% sm_fr_all={};
% for plotn = 1 : 4
%     for j=1:length(assigned_tastes)
%         fr_per_taste=[];
%         norm_fr_per_taste=[];
%           if norm==1
%                norm_fr_per_taste= FRMatrix{plotn,j};
%         fr = mean(norm_fr_per_taste,1);
%          elseif norm==0
%              fr_per_taste= FRMatrix{plotn,j};
%               fr = mean(fr_per_taste,1);
%          end
%         sm_fr= gaussmooth(fr,1000,2); % smoothing the firing rate
%         sm_fr_all{j,plotn}=sm_fr;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
fontSize = 12;
font = 'Arial';
timept=linspace(Pre/1000,Post/1000,(Post-Pre)/Binsize+1);
timepoint=timept(1:end-1)+Binsize/1000/2;
M=cell2mat(cellfun(@max,sm_fr_all,'uni',false));
Ymax=max(max(M,[],'omitnan'));
N=cell2mat(cellfun(@min,sm_fr_all,'uni',false));
Ymin=min(min(N,[],'omitnan'));
colors=distinguishable_colors(8);
if contains(foldername, 'sucrose')
    colors(6,:)=[];%generate some vectors contaning color values
else
end
h=gobjects(2,4);
t=tiledlayout(2,2);
subPlotNames={'increased firing Q331K', 'decreased firing Q331K' , 'increased firing non-Tg' ,'decreased firing non-Tg'};
% subPlotNames={['increased firing Q331K (n=',num2str(size(list{1,1},1)),')'],['decreased firing Q331K (n=',num2str(size(list{1,2},1)),')'] ,['increased firing non-Tg (n=',num2str(size(list{1,3},1)),')'] ,['decreased firing non-Tg (n=',num2str(size(list{1,4},1)),')']};
for plotn=1:4
    ax(plotn)=nexttile;
    plotsubtitle=[];
    for j = 1:length(assigned_tastes)
        if ~isnan(sm_fr_all{j,plotn})
            h(:,plotn)=plot(timepoint,sm_fr_all{j,plotn});hold on;
            set(h(:,plotn),'LineWidth',2,'Color',colors(j,:));
            %set(h(:,plotn),'LineWidth',2,'Color',color_gradient(j,:));
            hold on
            plotsubtitle=[plotsubtitle, num2str(size(FRMatrix{plotn,j},1)),' '];
        else
            plotsubtitle=[plotsubtitle, '0 '];
            continue
        end
    end
    xticks(ax(plotn),[xlim_start/1000:2: xlim_finish/1000]);
    title(subPlotNames{plotn}, plotsubtitle,'FontSize', fontSize,'FontName',font)
    plot([0 0],get(gca,'ylim'),'k--')
    box(ax(plotn),'off');pbaspect([1 1 1]);
    xtickangle(ax(plotn),0);
end
if norm==1
    ylim(ax(1),[0.45,Ymax]);
    ylim(ax(3),[0.45,Ymax]);
    ylim(ax(2),[0.25,0.55]);
    ylim(ax(4),[0.25,0.55]);
else
    ylim(ax(plotn),[Ymin,Ymax])
end
xlim(ax,[xlim_start/1000 xlim_finish/1000]);
xlabel(t,'Time(s)','FontSize', fontSize,'FontName',font);
if norm==1
    ylabel(t,'Normalized firing rate','FontSize', fontSize,'FontName',font)
elseif norm==0
    ylabel(t,'Firing rate (Hz)','FontSize', fontSize,'FontName',font)
end
legend(ax(2),[assigned_tastes,strrep(stimulus,'_',' ')],'location','northeastoutside','FontSize', fontSize,'FontName',font);
set(ax,'TickDir','out','color','none','box','off');
set(t,'units','normalized','outerposition',[0 0 1 1])
t.TileSpacing = 'compact';
t.Padding = 'compact';
% f=gcf;
% f.Position=get(0,'ScreenSize');
%% save files
if norm==1
    foldername= [foldername,'_double_removed_auROC'];
elseif norm==0
    foldername= [foldername,'_double_removed_nonorm'];
end
if baseline_option==2
    subfoldername=['bw_sd',num2str(bw_sd),'_baseline_pre_stimulus_',num2str(baseline_window),'_binsize_',num2str(Binsize),'ms_xlim_',num2str(-xlim_start),'_',num2str(xlim_finish),'_pre',num2str(-Pre),'_post',num2str(Post),'_prewindow',num2str(pre_window)];
elseif baseline_option==1
    subfoldername=['bw_sd',num2str(bw_sd),'_baseline_',num2str(baseline_window),'_binsize_',num2str(Binsize),'ms_xlim_',num2str(-xlim_start),'_',num2str(xlim_finish),'_pre',num2str(-Pre),'_post',num2str(Post),'_prewindow',num2str(pre_window)];
end
cd(frsavepath);
if exist(foldername,"dir")==7
    cd([frsavepath,'\',foldername]);
else
    mkdir(foldername);
end
cd([frsavepath,'\',foldername]);
if exist(subfoldername,"dir")==7
    cd([frsavepath,'\', foldername,'\',subfoldername]);
else
    mkdir(subfoldername);
end
cd([frsavepath,'\', foldername,'\',subfoldername]);
saveas(gcf,figurename);
exportgraphics(gcf,[figurename,'.pdf'],ContentType="vector");
save('FRMatrix_per_plot_taste','FRMatrix');
%%%%%%%%%%%%%%ANOVA for FRmatrix%%%%%%%%%%%%%%%
trialtype=assigned_tastes;
Timebin_Labels=num2cell((anova_start/Binsize+1:(anova_end-anova_start)/Binsize));
%for excitatory responses
temp=[];
for i=1:size(FRMatrix,1)
    temp=squeeze(mean(reshape(vertcat(FRMatrix{i,:}),size(assigned_tastes,2),[],(-Pre+Post)/Binsize),1));
    avgFRmatrix{i}=temp(:,(anova_start-Pre)/Binsize+1:(anova_end-Pre)/Binsize);
    nTrial=size(temp,1);
    Timebin_temp=repmat(Timebin_Labels,nTrial,1);
    Timebin_m{i}=reshape(Timebin_temp,[],1);
end
data=vertcat(avgFRmatrix{1}(:),avgFRmatrix{3}(:)); % avgFRmatrix: 1 x 4 cell (Q331K increase, Q decrease, nonTg imcrease, nonTg decrease) each with # of neurons x #of timebins double
times=cell2mat(vertcat(Timebin_m{1},Timebin_m{3}));
g_Q331K =repmat({'Q331K'},size(avgFRmatrix{1},1)*size(avgFRmatrix{1},2),1);
g_nonTg =repmat({'nonTg'},size(avgFRmatrix{3},1)*size(avgFRmatrix{3},2),1);
genotypes=char(vertcat(g_Q331K,g_nonTg));
[pval,tbl,stats] = anovan(data,{times genotypes},'Model','full','varnames',{'times', 'genotypes'},'display','on');
% Convert table to string
tblStr = formattedDisplayText(tbl);
% Write string to file
fid = fopen('anovaTable_excitatory.txt', 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')
save("ANOVA_averaged_across_taste_excitatory",'stats','pval');
% for inhibitory responses
data=vertcat(avgFRmatrix{2}(:),avgFRmatrix{4}(:));
times=cell2mat(vertcat(Timebin_m{2},Timebin_m{4}));
g_Q331K =repmat({'Q331K'},size(avgFRmatrix{2},1)*size(avgFRmatrix{2},2),1);
g_nonTg =repmat({'nonTg'},size(avgFRmatrix{4},1)*size(avgFRmatrix{4},2),1);
genotypes=char(vertcat(g_Q331K,g_nonTg));
[pval,tbl,stats] = anovan(data,{times genotypes},'Model','full','varnames',{'times', 'genotypes'},'display','on');
% Convert table to string
tblStr = formattedDisplayText(tbl);
% Write string to file
fid = fopen('anovaTable_inhibitory.txt', 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')
save("ANOVA_averaged_across_taste_inhibitory",'stats','pval');
end