% LFP script to view the spectrogram and correlation coefficient plot for
% one channel across different sessions at the same time to look for
% consistency of the channel spectrograms over dates

% Initialises session and channels
for channel=[6,11,17,21,26,30]
    
    sessionno="session01";
    arrayname="array01";
    sessionname = {'20181001', '20181004', '20181005', '20181008', '20181009',...
        '20181010','20181011', '20181015', '20181016', '20181017', '20181022', '20181026'};
    
    % Initialises channel name
    if numel(num2str(channel))==1
        channelname=strcat("channel00",string(channel),"");
    elseif numel(num2str(channel))==2
        channelname=strcat("channel0",string(channel),"");
    elseif numel(num2str(channel))==3
        channelname=strcat("channel",string(channel),"");
    end
    
    % Creates a new structure which stores the relevant data for the channel in
    % each different session
    for sessionname=sessionname
        % change directory to each day
        folder=strcat("/Volumes/Hippocampus/Data/picasso/",sessionname,"/",sessionno,"");
        cd(folder)
        % load the combined matrix for all channels
        load('arrayspec.mat');
        % retrieve the data for the required channel
        chn.(strcat("Date_",sessionname,""))=arrayspec(channel);
        
        % clear existing arrayspec structure
        clear arrayspec;
    end
    
    % -------------------------------------------------------------------------
    % plot the spectrogram for all sessions with the same channel
    
    figure("Position", get(0, "Screensize"))
    subplotno=1;
    
    for sessionname=sessionname
        % Initialises subplot
        subplot(4,3,subplotno);
        
        % Generates subplot for spectrogram
        surf(chn029.(strcat("Date_",sessionname,"")).T,chn029.(strcat("Date_",sessionname,"")).F,chn029.(strcat("Date_",sessionname,"")).Pnorm,"EdgeColor","none");
        axis xy; axis([-0.5 inf 0 150]); colormap(jet); view(0,90); caxis([-3 3]);
        title(strcat("session",string(sessionname),""),"FontSize",6);
        set(gca,"FontSize",6);
        
        % Plot lines to mark the cue presentation period
        hold on
        line([0 0],[0 150],[100 100; 100 100],'Color','k');
        line([1 1],[0 150],[100 100; 100 100],'Color','k');
        hold off
        
        subplotno=subplotno+1;
    end
    
    colorbar("Position", [0.915  0.11  0.015  0.8]);
    sgtitle(strcat("Normalised PSD averaged for all trials in each session for ",channelname,""),"FontSize",12);
    
    % Save the figure
    saveas(gcf,strcat("/Users/AlizarinMoon/Dropbox/MATLAB/Hippocampal LFP Research/Theta/",channelname,"_spec.png",""))
    close
    
    % -------------------------------------------------------------------------
    % plot the non-standardised cc plot for all sessions with the same channel
    
    figure("Position", get(0, "Screensize"))
    subplotno=1;
    
    for sessionname=sessionname
        % Initialises subplot
        subplot(4,3,subplotno);
        
        % Generates subplot for cc plot
        heatmap(chn029.(strcat("Date_",sessionname,"")).F,flipud(chn029.(strcat("Date_",sessionname,"")).F),flipud(chn029.(strcat("Date_",sessionname,"")).cc),...
            "XLimits",{0,148.43750},"YLimits",{148.43750,0},...
            "XDisplayLabels",(round(chn029.(strcat("Date_",sessionname,"")).F)),"YDisplayLabels",(flipud(round(chn029.(strcat("Date_",sessionname,"")).F))),...
            "Colormap", jet,"ColorLimits",[-0.2,0.7],"ColorbarVisible","on",...
            "Title",strcat("session",string(sessionname),""),"FontSize",5);
        subplotno=subplotno+1;
    end
    
    sgtitle(strcat("CC plot averaged for all trials in each session for ",channelname,""),"FontSize",12);
    
    % Save the figure
    saveas(gcf,strcat("/Users/AlizarinMoon/Dropbox/MATLAB/Hippocampal LFP Research/Theta/",channelname,"_cc.png",""))
    close
    
end

% -------------------------------------------------------------------------

% % initialise spectrogram subplot for all sessions with channel 29
% figure("Position", get(0, "Screensize"))
% subplotno=1;
% 
% for sessionname=sessionname
%     % change directory to each day
%     folder=strcat("/Volumes/Hippocampus/Data/picasso/",sessionname,"/",sessionno,"");
%     cd(folder)
%     % load the combined matrix for all channels
%     load('arrayspec.mat');
%     
%     % Initialises subplot
%     subplot(4,3,subplotno);
%     
%     % Generates subplot for spectrogram
%     surf(arrayspec(29).T,arrayspec(29).F,arrayspec(29).Pnorm,"EdgeColor","none");
%     axis xy; axis([-0.5 inf 0 150]); colormap(jet); view(0,90); caxis([-3 3]);
%     title(strcat("session",string(sessionname),""),"FontSize",6);
%     set(gca,"FontSize",6);
%     subplotno=subplotno+1;
%     
%     % clear existing arrayspec structure
%     clear arrayspec;
% end
% 
% colorbar("Position", [0.915  0.11  0.015  0.8]);
% sgtitle(strcat("Normalised PSD averaged for all trials in each session for ",channelname,""),"FontSize",12);
% 
% % Save the figure
% saveas(gcf,strcat("/Users/AlizarinMoon/Dropbox/MATLAB/Hippocampal LFP Research/Theta/",channelname,".png",""))
% close

