load("202210_4taste_celltype_summary.mat")
p=tiledlayout(1,4);
hg=gobjects(1,4);
hg(1,1)=nexttile;
name='sig009a_wf';
file='20220326_T67_4taste_sorted.nex';
cl=lines(8);
colors=cl(6,:); % blue, pyramidal
type=1; %type =1: print 1 window (out of 4) where the trough is
%type =2; %print 4 windows
% extract waveforms for .nex and print
namechunks = strsplit(file,'_');
date = namechunks{1};
mouseID = namechunks{2};
Nex= readNexFile(file);
% find the waveform index number
for i=1:size(Nex.waves,1)
    if string(Nex.waves{i,1}.name)==string(name)
        wf=Nex.waves{i,1}.waveforms;
    else
    end
end
filename=['waveform_',date,'_',mouseID,'_',name];
fontSize=12;
font='Arial';
plot_waveform_f(wf,filename,colors,type);  hold on;
name='sig017d_wf';
file='20220705_T88_4taste_sorted.nex';
cl=lines(8);
colors=cl(2,:); % red, interneuron
namechunks = strsplit(file,'_');
date = namechunks{1};
mouseID = namechunks{2};
Nex= readNexFile(file);
% find the waveform index number
for i=1:size(Nex.waves,1)
    if string(Nex.waves{i,1}.name)==string(name)
        wf=Nex.waves{i,1}.waveforms;
    else
    end
end
filename=['waveform_',date,'_',mouseID,'_',name];
fontSize=12;
font='Arial';
plot_waveform_f(wf,filename,colors,type);
ylim([-0.08 0.06]);
yticks(-0.08:0.02:0.06);
pbaspect([1 1 1])
axis square
ylabel('Amplitude (µV)')
xlabel('Time (µs)')
hg(1,2)=nexttile;
histogram(C_new.Trough2Peak(C_new.T2Pcelltype==1),'FaceAlpha',0.8,'FaceColor',cl(6,:),'EdgeColor', cl(6,:), 'BinWidth',1e-2) %excitatory non_Tg
hold on
histogram(C_new.Trough2Peak(C_new.T2Pcelltype==2),'FaceAlpha',0.8,'FaceColor',cl(2,:),'EdgeColor', cl(2,:),'BinWidth',1e-2) %excitatory non_Tg
legend({'putative pyrimidal','putative interneuron'},"Box","off",Location="northwest")
axis square
ylabel('# of units')
xlabel('Trough-to-peak (µs)')

hg(1,3)=nexttile;
fontSize=12;
font='Arial';
% calculate average baseline firing rates based on T2Pcelltype
avg_baseline_T2Pcelltype=zeros(2,max(C_new.T2Pcelltype));
ste_baseline_T2Pcelltype=zeros(2,max(C_new.T2Pcelltype));
baseline_T2Pcelltype={};
for k=1:max(C_new.T2Pcelltype)
    i1=(C_new.mutant==1&C_new.T2Pcelltype==k);
    i2=(C_new.mutant==0&C_new.T2Pcelltype==k);
    avg_baseline_T2Pcelltype(1,k)=mean(C_new.avg_baseline(i1),1);
    avg_baseline_T2Pcelltype(2,k)=mean(C_new.avg_baseline(i2),1);
    ste_baseline_T2Pcelltype(1,k)=std(C_new.avg_baseline(i1),[],1)/sqrt(sum(i1));
    ste_baseline_T2Pcelltype(2,k)=std(C_new.avg_baseline(i2),[],1)/sqrt(sum(i2));;
    baseline_T2Pcelltype{1,k}=C_new.avg_baseline(i1); %Q331K
    baseline_T2Pcelltype{2,k}=C_new.avg_baseline(i2); %nonTg
end
% cumulative distribution for baseline firing rate (non_Tg first)
cl=distinguishable_colors(8);
reshaped_baseline= reshape(baseline_T2Pcelltype,4,1); %Q331K Pyr, nonTg Pyr, Q331K int, nonTg int
colors_reshaped=vertcat(cl(5,:),cl(4,:),cl(5,:),cl(4,:));
f=cell(1,4);x=cell(1,4);
for i=1:size(reshaped_baseline,1)
    [f{i},x{i}]=ecdf( reshaped_baseline{i,1});
    h(i)=plot(x{i},f{i},'color',colors_reshaped(i,:),'LineStyle','-','LineWidth',2);
    if i>2
        set(h(i),'LineStyle','--')
    end
    hold on
end
legend ({'Q331K Putative Pyramidal','non-Tg Putative Pyramidal','Q331K Putative Interneuron','non-Tg Putative Interneuron'},'color','none','box','off')
ylabel('Frequency');
xlabel('Firing rate (Hz)')
axis square

hg(1,4)=nexttile;
% xvals = bar_custom(avg_licks_by_indv_trials');
h0=bar(avg_baseline_T2Pcelltype','FaceColor','flat');hold on
colors=distinguishable_colors(6);
set(h0(:,1),'FaceColor',colors(5,:),'EdgeColor',colors(5,:));
set(h0(:,2),'FaceColor',colors(4,:),'EdgeColor',colors(4,:));
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg_baseline_T2Pcelltype');
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg_baseline_T2Pcelltype(i,:), ste_baseline_T2Pcelltype(i,:), 'k', 'linestyle', 'none');
end
hold off
xticklabels({'Putative Pyr.','Putative Inh.'})
ylabel('Firing rate (Hz)')
axis square 
set(hg,'TickDir','out','box','off');
set(hg,'FontSize', fontSize,'FontName',font);
f=gcf;
f.Position=get(0,'ScreenSize');
figurename=['Figure2.pdf'];
exportgraphics(gcf,figurename)