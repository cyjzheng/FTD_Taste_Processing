function [sig_anova]=fun_anova_datadf(datadf_PSTH, assigned_tastes, Binsize, start_window, end_window, Pre, alpha,path,varargin)
% alpha=0.05;
% path='D:\AF Lab\matlab_data\best_window\202210_4taste_pre6000post10000goodcentral_v2b1\SD2\bw_avgALL_responsive_valve_on_pre1000_post250_postwindow250_limit6000';
% assigned_tastes={'S','N','Q','C'};
% Binsize=250;
% Pre=-6000;
% end_window=6000;
% start_window=0;
%varargin{1} is meet_fr_threshold; if entered, will adjust the # total units for proportion_final
%varargin{end} is sig_anova
trialtype=assigned_tastes;
Timebin_Labels=num2cell((start_window-Pre)/Binsize+1:(end_window-Pre)/Binsize);
nNeuron=size(unique(datadf_PSTH.neuron_numeric_ID),1);
% if size(varargin,2)<1||size(varargin,2)<2&&isa(varargin{end},'logical')
    disp('calculating sig_anova')
for n=1:nNeuron
    n1_df=datadf_PSTH(datadf_PSTH.neuron_numeric_ID==n&datadf_PSTH.Taste_numeric_ID<=size(assigned_tastes,2),:);
    all_scmatrix=[cell2mat(n1_df.(append('scmatrix_',num2str(Binsize))))];
    temp=all_scmatrix(:,(start_window-Pre)/Binsize+1:(end_window-Pre)/Binsize);
    data=reshape(temp,[],1);
    nTrial=size(all_scmatrix,1);
    %create labels for time bins
    %make a matrix with similar dimension as the spike count matrices
    Timebin_temp=repmat(Timebin_Labels,nTrial,1);
    times=cell2mat(reshape(Timebin_temp,[],1));
%     taste_labels =cell2mat( [n1_df.assigned_tastes{:}]');
    taste_labels = string(n1_df.assigned_tastes);
    taste = repmat(taste_labels,(-start_window+end_window)/Binsize,1);
    [pval{n},tbl{n},stats{n}] = anovan(data,{times taste},'Model','full','varnames',{'times', 'taste'},'display','off');
end
temp_p=horzcat(pval{1,:});
t=(temp_p(2,:)<alpha|temp_p(3,:)<alpha)';%either the taste(trial type) main effect or the taste x times interaction is significant
sig_anova=find(t)';
% elseif ~isa(varargin{end},'logical')
%   disp('sig_anova found')
%     sig_anova=varargin{end};
% end 
cd(path)
f1=load('sig_bw');
t1=f1.sig_bw;
[v,~]=intersect(t1(:,1),sig_anova);
%then trying to find the matched rows with repeats in sig_bw(t1)
b_sig = ismember(t1(:,1),v);
i_sig = find(b_sig);
anova_sig_bw=t1(i_sig,:);
%find unique
f2=load('uni_sig_bw');
t2=f2.uni_sig_bw;
[~,i_kw_t2]=intersect(t2(:,1),sig_anova);
anova_uni_sig_bw=t2(i_kw_t2,:);
%% find responsive per taste using sig
resp_by_taste=anova_sig_bw;
%final: (# tastes x increase (1)/ decrease (2)x mutant (mutant:1,
%control:2) )
anova_final=zeros(max(resp_by_taste(:,2)),2,2);
%    for j=1: length(tastes) % for each taste
for j=1: max(resp_by_taste(:,2))
    temp1=resp_by_taste(:,2)==j;
    temp2= resp_by_taste(temp1(:),:);
    idx=(temp2(:,3)==1 & temp2(:,4)==1); % increase (column 3) and mutant (column 4)
    anova_final(j,1,1)=sum(idx(:));
    idx2=(temp2(:,3)==0 & temp2(:,4)==1); % decrease (column 3) and mutant (column 4)
    anova_final(j,2,1)=sum(idx2(:));
    idx3=(temp2(:,3)==1 & temp2(:,4)==0); % increase (column 3) and control (column 4)
    anova_final(j,1,2)=sum(idx3(:));
    idx4=(temp2(:,3)==0 & temp2(:,4)==0); % decrease (column 3) and control (column 4)
    anova_final(j,2,2)=sum(idx4(:));
end
[~,id]=unique(datadf_PSTH.neuron_numeric_ID);
m=table2array(datadf_PSTH(id,'mutant'));
if size(varargin,2)<1
       disp('no meet_fr_threshold found')
n_Q331K=nnz(m);
n_nonTg=(size(m,1)-nnz(m));
else  
   disp('meet_fr_threshold found')
    meet_fr_threshold= varargin{1,1};
n_Q331K=nnz(m==1&meet_fr_threshold');
n_nonTg=nnz(m==0&meet_fr_threshold');
end
anova_proportion_final(:,:,1)=anova_final(:,:,1)/n_Q331K;
anova_proportion_final(:,:,2)=anova_final(:,:,2)./n_nonTg;
folder=['anova_bin',num2str(Binsize)];
if exist(folder, 'dir')
    cd([path,'\' folder]);
else
    mkdir(folder);
    cd([path,'\' folder]);
end
save('anova_sig','sig_anova')
save('anova_sig_bw',"anova_sig_bw");
save(['anova_uni_sig_bw'],'anova_uni_sig_bw')
save(['anova_proportion_final'],'anova_proportion_final');
save(['anova_final'],'anova_final');

end