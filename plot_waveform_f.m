function [p, h]=plot_waveform_f(wf,name,colors,type)
% inputs:
% name: the name of the waveform, a string
% wf: a 128 x n double(n varies depending on recording length; 128 is the y dimension for all waveforms) exported from .nex file
av_w=mean(wf,2)';
std_dev2=std(wf,0,2)';
Xlim=1:length(av_w);
%making a shadow with standard error
xx = [Xlim(1:end), fliplr(Xlim(1:end))];
% yy = [av_w + std_dev2/sqrt(size(wf,1)),...
%     fliplr(av_w - std_dev2/sqrt(size(wf,1)))];

%just plot the st dev
yy = [av_w + std_dev2,...
    fliplr(av_w - std_dev2)];
h = fill(xx, yy, colors,'LineStyle','none');
hold on 
p=plot(av_w,'w','LineWidth',1);
box('off');
% in order to plot the waveform with the largest amplitude, figure out
% which channel the trough is in
 [~,trough_ind] = min(av_w); % Get trough and its index
    n=0;
while n * 32 <trough_ind
    n=n+1;
end   
switch type
    %print 1 window (out of 4) where the trough is 
case 1
xlim([32*(n-1) 32*n]);
xticks(32*(n-1):4:32*n);
xticklabels(0:4*25:32*25);
pbaspect([0.4 1 1]) 
 %  print full waveforms (4 tetrode windows)
case 2
 xlim([0 128]);
pbaspect([2 1 1]) 
end

end

