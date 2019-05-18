function [S,cc] = lfpfig(spec, varargin)
% LFPFIG Displays the normalised PSD and correlation coefficient data
% contained within a spec structure as figures
%
%   Inputs: lfpfig(spec,type,number)
%   (1) spec = spectrogram stucture where spec.Pnorm contains the normalised
%   PSD data and spec.cc contains the correlation coefficient data
%   (2) type = "spec" or "mspec", indicating structure consists of individual
%   trials, or one structure averaged across all trials
%   (3) number = [] array consisting of the individual trial
%   numbers to plot if spec is a structure consisting of individual trials
%
%   Outputs:
%   S is the spectrogram handle
%   cc is the corr coeff plot handle

% =========================================================================

% checks for the input arguments
narginchk(1,3);
type=varargin{1};

% =========================================================================

% If structure type is the individual trials
if type=="spec"
    number=varargin{2};
    % Plots the normalized PSD data in a figure
    figure('Position', get(0, 'Screensize'))
    
    plot=1;
    for n=number
        subplot(3,3,plot);
        S=surf(spec(n).T,spec(n).F,spec(n).Pnorm,'EdgeColor','none');
        axis xy; axis([-0.5 inf 0 150]); colormap(jet); view(0,90); caxis([-10 10]);
        set(gca,'FontSize',6); xticks(-0.5:0.5:spec(n).T(end)); yticks(0:10:150);
        title(strcat("Normalised Spectrogram of Trial:",string(n)));
        plot=plot+1;
    end
    
    % Plot lines to mark the cue presentation period
    hold on
    line([0 0],[0 150],[10 10; 10 10],'Color','k');
    line([1 1],[0 150],[10 10; 10 10],'Color','k');
    hold off
    
    % colorbar for spectrogram plot
    colorbar;
    cbpos=get(subplot(3,3,9),'Position');
    colorbar('Position', [cbpos(1)+cbpos(3)+0.01  cbpos(2)  0.01  cbpos(4)*3.7]);
    
    % Plots 9 trials worth of CC plots in one figure
    figure('Position', get(0, 'Screensize'))
    
    plot=1;
    for n=number
        subplot(3,3,plot);
        cc=heatmap(spec(n).F,flipud(spec(n).F),flipud(spec(n).cc),...
            'XLimits',{0,150.39063},'YLimits',{150.39063,0},...
            'XDisplayLabels',(round(spec(n).F)),'YDisplayLabels',(flipud(round(spec(n).F))),...
            'XLabel','Frequency (Hz)','YLabel','Frequency (Hz)','FontSize',6,...
            'Colormap', jet,'ColorLimits',[-0.2,0.7],'ColorbarVisible','on',...
            'Title',strcat("Correlation coefficients of Trial:",string(n)));
        plot=plot+1;
    end
end         

% =========================================================================

% If structure type is the averaged trials
if type=="mspec"
    
    figure('Position', get(0, 'Screensize'))
    
    subplot(2,1,1); % Plots the normalized PSD data
    surf(spec.T,spec.F,spec.Pnorm,'EdgeColor','none');
    axis xy; axis([-0.5 inf 0 150]); colormap(jet); view(0,90); caxis([-3 3]);
    set(gca,'FontSize',6); xticks(-0.5:0.5:7); yticks(0:10:150);
    title("Normalised PSD averaged for all trials",'FontSize',10);
    colorbar;
    
    % Plot lines to mark the cue presentation period
    hold on
    line([0 0],[0 150],[100 100; 100 100],'Color','k');
    line([1 1],[0 150],[100 100; 100 100],'Color','k');
    hold off
    
    subplot(2,2,3); % Plots the CC plot
    heatmap(spec.F,flipud(spec.F),flipud(spec.cc),...
        'XLimits',{0,150.39063},'YLimits',{150.39063,0},...
        'XDisplayLabels',(round(spec.F)),'YDisplayLabels',(flipud(round(spec.F))),...
        'XLabel','Frequency (Hz)','YLabel','Frequency (Hz)','FontSize',6,...
        'Colormap', jet,'ColorLimits',[-0.2,0.7],'ColorbarVisible','on',...
        'Title','Correlation coefficient plot averaged for all trials');
end


