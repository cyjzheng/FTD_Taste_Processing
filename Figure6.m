%%  plot decoding accuracy
Pre=-6;
Post=10;
binsize=250;
filename=cell(1,2);
filename{1,1}='250ms_nonTg';
filename{1,2}='250ms_mutant';
colors=distinguishable_colors(8);
h=gobjects(1,2);
t=tiledlayout(1,2);
for i=1:2
    ax(i)=nexttile;
    hold on
    color_to_plot= colors(3+i,:);
    load([filename{i},'_basic_decoding_results_original_test_only_at_training_times']);
    %get mean decoding results for 20 runs
    temp=[DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS];
    temp1=horzcat(temp.mean_decoding_results);
    temp2=mean(temp1,2)';
    %get standard deviation
    curr_stdev=[temp.stdev];
    stdev=horzcat(curr_stdev.over_resamples);
    stdev=mean(stdev, 2)';
    Xlim=1:size(temp2,2);
    %making a shadow with standard deviation
    xx = [Xlim(1:end), fliplr(Xlim(1:end))];
    yy = [temp2 + stdev,...
        fliplr(temp2 - stdev)];
    colors=distinguishable_colors(8);
    p(:,1) = fill(xx, yy, color_to_plot,'LineStyle','none');
    set(p(:,1),{'DisplayName'},{'Original'});
    hold on
    p(:,2)=plot(temp2,'w','LineWidth',1);
    set(p(:,2),{'DisplayName'},{''});
    load([filename{i},'_basic_decoding_results_shuffled_test_only_at_training_times']);
    %get mean decoding results for 20 runs
    temp=[DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS];
    temp1=horzcat(temp.mean_decoding_results);
    temp2=mean(temp1,2)';
    %get standard deviation
    curr_stdev=[temp.stdev];
    stdev=horzcat(curr_stdev.over_resamples);
    stdev=mean(stdev, 2)';
    Xlim=1:size(temp2,2);
    %making a shadow with standard deviation
    xx = [Xlim(1:end), fliplr(Xlim(1:end))];
    yy = [temp2 + stdev,...
        fliplr(temp2 - stdev)];
    colors=distinguishable_colors(8);
    p(:,3)= fill(xx, yy, color_to_plot,'LineStyle','none','FaceAlpha',0.5);
    set(p(:,3),{'DisplayName'},{'Shuffled'});
    hold on
    p(:,4)=plot(temp2,'w','LineWidth',1);
    box('off');
    set(p(:,4),{'DisplayName'},{''});
    legend(ax(i),'location','northeastoutside');
    pbaspect([1 1 1])
    % axis square
end
linkaxes(ax,'x');
% set y axis
linkaxes(ax,'y');


xlim(ax,[(-2-Pre)*1000/binsize (5-Pre)*1000/binsize]);
xticks(ax,linspace((-2-Pre)*1000/binsize,(5-Pre)*1000/binsize,8));
xticklabels(ax,num2str(linspace(-2,5,8)'))
xlabel(t,'Time (s)','FontSize',12)
ylabel(t,'Decoding Accuracy','FontSize',12)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(ax,'TickDir','out','box','off');
cd(data_directory_name)
exportgraphics(t,[filename{1},'_',filename{2},'_','decoding_accuracy_plot_xlim.pdf'],'Resolution', 300)
