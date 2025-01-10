function [sig_bw,uni_sig_bw,inc,dec,latency]=best_window_final(assigned_tastes,stimulus,window,pre_window,post_window,Pre,step, limit,alpha,n_sd,datadf_PSTH,datasavepath,subfoldername)
%input: data_dfPSTH
% applies a moving window on the firing rate averaged across all trial types
%this version compares the response during the window (y_fr) to the
%baseline of that neuron for ALL trialtypes
%different from v2b1: more stringent- only those significant for 8 steps
%are selected
% SD is performed on the smoothed firing rate in pre_window, pre_window starts at -4s
% subfoldername='202210_4taste_pre6000post10000goodcentral_v2b1a';
% assigned_tastes={'S','N','Q','C'};
% stimulus='valve_on';
% window=250;
% pre_window=1000;
% post_window=250;
% Pre=-6000;
% Post=10000;
% step=25;
% limit=5000;
% alpha=0.05;
% n_sd=2; % set criterion for how good the best window is. response is n_sd*standard deviation compared to the baseline
alpha=alpha/length(assigned_tastes);
filename=['bw_avgALL_responsive_',stimulus,'_pre',num2str(pre_window),'_post',num2str(window),'_postwindow',num2str(post_window),'_limit',num2str(limit)];%bw=best window
%% make a moving sum for time 0 to limit, and compare that to the average baseline and SD for that trial type
neuronIDs=unique(datadf_PSTH.neuron_numeric_ID);
nNeurons=length(neuronIDs);
for i= 1:nNeurons
    n_df=datadf_PSTH(datadf_PSTH.neuron_numeric_ID==neuronIDs(i),:);
    scmatrix=cell2mat(n_df.(append('scmatrix_',num2str(step))));
    x_scmatrix=scmatrix(:,(-Pre-4000)/step+1:(-Pre-4000+pre_window)/step); % baseline starts from -4s before taste delivery
    baseline_fr= x_scmatrix/(step/1000);
    sm_baseline_fr=gaussmooth(mean(baseline_fr,1),2,100);
    avg_baseline(i)=mean(sm_baseline_fr,'all');
    SD(i)=std(sm_baseline_fr, 0, "all");
    for j=1:length(assigned_tastes)
        n1_df=datadf_PSTH(datadf_PSTH.neuron_numeric_ID==neuronIDs(i)&strcmp(string(datadf_PSTH.assigned_tastes),assigned_tastes{j}),:);
        movewindow=[];
        nTrials=size(n1_df,1);
        nWindows=1+(limit-window)/step;
        sc_matrix_best=[];
        %movewindow: a nTrials x nWindows matrix. each element is the spike
        %count for specific window during a trial for that trial type
        all_scmatrix=cell2mat(n1_df.(append('scmatrix_',num2str(step))));
        window_raw=  all_scmatrix(:,(-Pre/step+1):(-Pre+limit)/step);
        movewindow=movsum(window_raw,window/step,2,'Endpoints','discard');
        %         smoothed_window=gaussmooth(mean_window_raw,100,2);
        %         smooth_movewindow=movsum(smoothed_window,window/step,2,'Endpoints','discard');
        Fr_mean_resp=sum(movewindow,1)/nTrials/(window/1000); % calculate the firing rate during that window
        %substract the baseline from the average, check if the difference is larger than n_sd*SD
        FR_mean_resp{i,j}=Fr_mean_resp;
        Fr_diff=Fr_mean_resp-avg_baseline(i);
        FR_diff{i,j}=Fr_diff;
        % i_1: an array of indices for windows that meet the criteria
        [~,i1]=find(abs(Fr_diff)>n_sd*SD(i));
        [T_on ,T_off]=Timing_onset_offset(abs(Fr_diff)>n_sd*SD(i), 1:1:size(Fr_diff,2), 0.8,1,0);
        %find the max firing rate in all these windows that meet the criteria
        [~,i2]=max(abs(Fr_diff(abs(Fr_diff)>n_sd*SD(i))));
        %         if ~isempty(i_1) && ~isempty(i_2) && ismember(i_1(i_2)-1,i_1)==1 &&  ismember(i_1(i_2)+1,i_1)==1
        numbins_post_window=post_window/step;
        if ~isempty(i1) && ~isempty(i2)&& any(T_off-T_on>=8)
            if  i1(i2)+numbins_post_window+(-Pre/step-1) <(-Pre+limit)/step
                y=all_scmatrix(:,i1(i2)+(-Pre/step):i1(i2)+numbins_post_window+(-Pre/step)-1); % a k X post_window/step double of firing rates (k: # of trials for that trial type)
            elseif  i1(i2)+numbins_post_window+(-Pre/step-1) >(-Pre+limit)/step
                y=all_scmatrix(:,i1(i2)+(-Pre/step):(-Pre+limit)/step); % a k X post_window/step double of firing rates (k: # of trials for that trial type)
            end
            latency(i,j)=(i1(i2)-1)*step+0.5*window;
            [p0,~,~]=ranksum(x_scmatrix(:), y(:));
            pv(i,j)=p0;
            inc(i,j)=double(mean(y,'all')>mean(x_scmatrix,'all'));
            dec(i,j)=double(mean(y,'all')<mean(x_scmatrix,'all'));
            max_response(i,j)=Fr_diff(:,i1(i2));
            raw_max_response(i,j)=Fr_mean_resp(:,i1(i2));
            max_idx(i,j)=i1(i2);
        else
            sc_matrix_best(i,j,:)=NaN;
            pv(i,j)=NaN;
            inc(i,j)=NaN;
            dec(i,j)=NaN;
            latency(i,j)=NaN;
            max_response(i,j)=NaN;
            raw_max_response(i,j)=NaN;
            max_idx(i,j)=NaN;
        end
    end
end
%%
%change values in inc and dec to NaN if there is no sigificant results for
%those tastes 6/23/22
t =pv<alpha;
inc=double(inc);
inc(~t) = NaN;
dec=double(dec);
dec(~t) = NaN;
ind=find(pv<alpha); %a linear indices of signifiance values in pv
[r,c]=ind2sub(size(pv), ind); %convert linear indices to row and column indices
sig_bw=[];
sig_bw(:,:)=[r,c];
% column 3: increase for that taste Yes=1, No=0
for i=1:size(sig_bw,1)
    sig_bw(i,3)=inc(sig_bw(i,1),sig_bw(i,2));
end
%column 4: mutant Yes=1
[~,id]=unique(datadf_PSTH.neuron_numeric_ID);
m=table2array(datadf_PSTH(id,'mutant'));
sig_bw(:,4)=m(sig_bw(:,1),:);
%% find responsive per taste using sig
resp_by_taste=sig_bw;
%final: 7 x2 x2 (# tastes x increase (1)/ decrease (2)x mutant (mutant:1,
%control:2) )
final=zeros(max(resp_by_taste(:,2)),2,2);
%    for j=1: length(tastes) % for each taste
for j=1: max(resp_by_taste(:,2))
    temp1=resp_by_taste(:,2)==j;
    temp2= resp_by_taste(temp1(:),:);
    idx=(temp2(:,3)==1 & temp2(:,4)==1); % increase (column 3) and mutant (column 4)
    final(j,1,1)=sum(idx(:));
    idx2=(temp2(:,3)==0 & temp2(:,4)==1); % decrease (column 3) and mutant (column 4)
    final(j,2,1)=sum(idx2(:));
    idx3=(temp2(:,3)==1 & temp2(:,4)==0); % increase (column 3) and control (column 4)
    final(j,1,2)=sum(idx3(:));
    idx4=(temp2(:,3)==0 & temp2(:,4)==0); % decrease (column 3) and control (column 4)
    final(j,2,2)=sum(idx4(:));
end
proportion_final(:,:,1)=final(:,:,1)/nnz(m);
proportion_final(:,:,2)=final(:,:,2)./(size(m,1)-nnz(m));
%% max response
%max response for increased responese
inc_max=max_response.*inc;
inc_max(inc==0)=NaN;
%max response for decreased responese
dec_max=max_response.*dec;
dec_max(dec==0)=NaN;
% cd('C:\Users\camel\PycharmProjects\pythonProject\data\latency\');
% if exist(subfoldername)==7
%     cd(['C:\Users\camel\PycharmProjects\pythonProject\data\latency\' subfoldername]);
% else
%     mkdir(subfoldername);
%     cd(['C:\Users\camel\PycharmProjects\pythonProject\data\latency\' subfoldername]);
% end
% if exist(['SD',num2str(n_sd)])==7
%     cd(['C:\Users\camel\PycharmProjects\pythonProject\data\latency\' subfoldername, '\','SD',num2str(n_sd)])
% else
%     mkdir(['SD',num2str(n_sd)])
%     cd(['C:\Users\camel\PycharmProjects\pythonProject\data\latency\' subfoldername, '\','SD',num2str(n_sd)])
% end
% if exist(filename)==7
%     cd(['C:\Users\camel\PycharmProjects\pythonProject\data\latency\' subfoldername, '\','SD',num2str(n_sd),'\', filename])
% else
%     mkdir(filename)
%     cd(['C:\Users\camel\PycharmProjects\pythonProject\data\latency\' subfoldername, '\','SD',num2str(n_sd),'\', filename])
% end
% writematrix(y,'latency_by_taste.csv');
% writematrix(err,'err_by_taste.csv');
% writematrix(latency,'raw_latency.csv');
% writematrix(m,'mutant_info.csv');
% writecell(Labels,'assigned_tastes.txt');
%% some more stuff about latency
counts=sum(~isnan(latency),2);
%count the number of units that have a latency value
counts_by_genotype(1,:)=size(nonzeros(counts(m==1,:)),1);
counts_by_genotype(2,:)=size(nonzeros(counts(m==0,:)),1);
%%
cd(datasavepath);
if exist(subfoldername)==7
    cd([datasavepath,'\', subfoldername]);
else
    mkdir(subfoldername);
    cd([datasavepath,'\', subfoldername]);
end
if exist(['SD',num2str(n_sd)])==7
    cd([datasavepath,'\', subfoldername, '\','SD',num2str(n_sd)])
else
    mkdir(['SD',num2str(n_sd)])
    cd([datasavepath,'\', subfoldername, '\','SD',num2str(n_sd)])
end
if exist(filename)==7
    cd([datasavepath,'\', subfoldername, '\','SD',num2str(n_sd),'\', filename])
else
    mkdir(filename)
    cd([datasavepath,'\', subfoldername, '\','SD',num2str(n_sd),'\', filename])
end
save('sig_bw','sig_bw')
%unique responsive
% column 1: # in summaryData/PSTH
% column 2: increased for one or more tastes Yes=1, No=0
% column 3: decreased for one or more tastes Yes=1, No=0
[~,ia,~]=unique(sig_bw(:,1),'row');
uni_sig_bw=sig_bw(ia,:);
uni_sig_bw(:,2)=sum(inc(sig_bw(ia),:),2,'omitnan')>=1;
uni_sig_bw(:,3)=sum(dec(sig_bw(ia),:),2,'omitnan')>=1;
save('avg_baseline','avg_baseline')
save("SD","SD")
save("FR_diff","FR_diff")
save('pv','pv');
save('uni_sig_bw','uni_sig_bw');
save('dec','dec');
save('inc','inc');
save('raw_latency_by_taste','latency');
save('max_response','max_response');
save('raw_max_response','raw_max_response');
save('final','final');
save('proportional_final','proportion_final');
save('mutant_info','m');
save('counts_by_genotype','counts_by_genotype')
%% latency
%latency for increased responese
inc_latency=latency.*inc;
inc_latency(inc==0)=NaN;
%latency for decreased responese
dec_latency=latency.*dec;
dec_latency(dec==0)=NaN;
[h(1,:),p(1,:)]=plot_latency(latency,m,'all',assigned_tastes);
[h(2,:),p(2,:)]=plot_latency(inc_latency,m,'inc',assigned_tastes);
[h(3,:),p(3,:)]=plot_latency(dec_latency,m,'dec',assigned_tastes);
save('kstest_latency_per_taste','p','h');
end