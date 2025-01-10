function plot_firing_rate_two_genotypes (Pre, Post, Binsize, Q331K_fr,Q331K_ste, nonTg_fr, nonTg_ste, fontSize,font, xlim_start,xlim_finish, assigned_tastes,std_option)
timept=linspace(Pre/1000,Post/1000,(Post-Pre)/Binsize+1);
timepoint=timept(1:end-1)+Binsize/1000/2;
% a=[vertcat(sm_fr_all{:,1}),vertcat(sm_fr_all{:,3})];
% max(a,[],'all')
Ymax=0.05+max(Q331K_fr(:,(xlim_start-Pre)/Binsize+1:(xlim_finish-Pre)/Binsize) + Q331K_ste(:,(xlim_start-Pre)/Binsize+1:(xlim_finish-Pre)/Binsize),[],'all');
Ymin=min(Q331K_fr(:,(xlim_start-Pre)/Binsize+1:(xlim_finish-Pre)/Binsize) - Q331K_ste(:,(xlim_start-Pre)/Binsize+1:(xlim_finish-Pre)/Binsize),[],'all')-0.05;
colors=distinguishable_colors(8);
if size(Q331K_fr,1)==1
    if std_option==1
        %making a shadow with standard error
        figure;
        xx = [timepoint(1:end), fliplr(timepoint(1:end))];
        yy = [ Q331K_fr(1,:) +  Q331K_ste(1,:),fliplr( Q331K_fr(1,:) -  Q331K_ste(1,:))];
        % h = fill(xx, yy, colors(5,:),'LineStyle','none',FaceAlpha=0.5); hold on;
        h = fill(xx, yy, colors(5,:),'LineStyle','none'); hold on;
        p=plot(timepoint,Q331K_fr(1,:),'w','LineWidth',1);
        hold on
        yy2 = [nonTg_fr(1,:) + nonTg_ste(1,:),fliplr(nonTg_fr(1,:) - nonTg_ste(1,:))];
        %  h2 = fill(xx, yy2, colors(4,:),'LineStyle','none',FaceAlpha=0.5); hold on;
        h2 = fill(xx, yy2, colors(4,:),'LineStyle','none'); hold on;
        p2=plot(timepoint,nonTg_fr(1,:) ,'w','LineWidth',1);hold on
    else
        p=plot(timepoint,Q331K_fr(1,:),'w','LineWidth',3,Color=colors(5,:));
        hold on
        p2=plot(timepoint,nonTg_fr(1,:) ,'w','LineWidth',3,Color=colors(4,:));
    end
    set(gca,'FontSize', fontSize,'FontName',font);
    box('off');
    daspect([1 (Ymax-Ymin)/((xlim_finish-xlim_start)/1000) 1]);
    pbaspect([1 1 1]);
    %     xticks(ax(j),[xlim_start/1000:1: (xlim_finish+500)/1000]);
    if Ymin>0.4
        ylim([0.45,Ymax]); % if plotting increased firing
    else
        ylim([0.25,0.55]); % if plotting decreased firing
    end
    ylabel('Normalized firing rate','FontSize', fontSize,'FontName',font);
    xlim([xlim_start/1000 (xlim_finish)/1000]);
    xlabel('Time(s)','FontSize', fontSize,'FontName',font)

    set(gca,'TickDir','out','color','none','box','off');

elseif size(Q331K_fr,1)==size(assigned_tastes,2)
    t=tiledlayout(1,size(assigned_tastes,2));
    for j=1:size(assigned_tastes,2)
        ax(j)=nexttile;
        if std_option==1
            %making a shadow with standard error
            xx = [timepoint(1:end), fliplr(timepoint(1:end))];
            yy = [ Q331K_fr(j,:) +  Q331K_ste(j,:),fliplr( Q331K_fr(j,:) -  Q331K_ste(j,:))];
            % h = fill(xx, yy, colors(5,:),'LineStyle','none',FaceAlpha=0.5); hold on;
            h = fill(xx, yy, colors(5,:),'LineStyle','none'); hold on;
            p=plot(timepoint,Q331K_fr(j,:),'w','LineWidth',1);
            hold on
        else
            p=plot(timepoint,Q331K_fr(j,:),'w','LineWidth',3,Color=colors(5,:));
            hold on
        end
        if std_option==1
            yy2 = [nonTg_fr(j,:) + nonTg_ste(j,:),fliplr(nonTg_fr(j,:) - nonTg_ste(j,:))];
            %  h2 = fill(xx, yy2, colors(4,:),'LineStyle','none',FaceAlpha=0.5); hold on;
            h2 = fill(xx, yy2, colors(4,:),'LineStyle','none'); hold on;
            p2=plot(timepoint,nonTg_fr(j,:) ,'w','LineWidth',1);hold on
        else
            p2=plot(timepoint,nonTg_fr(j,:) ,'w','LineWidth',3,Color=colors(4,:));
        end
        axis equal;
        set(gca,'FontSize', fontSize,'FontName',font);
        box('off');
        daspect(ax(j),[1 (Ymax-Ymin)/((xlim_finish-xlim_start)/1000) 1]);
        pbaspect(ax(j),[1 1 1]);
        %     xticks(ax(j),[xlim_start/1000:1: (xlim_finish+500)/1000]);
        if Ymin>0.4
            ylim(ax(j),[0.45,Ymax]); % if plotting increased firing
        else
            ylim(ax(j),[0.25,0.55]); % if plotting decreased firing
        end
    end
    ylabel(t,'Normalized firing rate','FontSize', fontSize,'FontName',font);
    xlim(ax,[xlim_start/1000 (xlim_finish)/1000]);
    xlabel(t,'Time(s)','FontSize', fontSize,'FontName',font)
    % set x axis
    linkaxes(ax,'x');
    % set y axis
    linkaxes(ax,'y');
    set(ax,'TickDir','out','color','none','box','off');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end
end
