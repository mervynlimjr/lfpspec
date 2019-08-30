% LFP script to view the spectrogram and correlation coefficient plot for
% one channel across different sessions at the same time to look for
% consistency of the channel spectrograms over dates

% Dropbox user folder
userfol="AlizarinMoon";

% Initialises session, channels, subplot variables
sessionno="session01";
arrayname="array01";
sessiondate=["20180221","20180312","20180323","20180402","20180404","20180416","20180625","20180704","20180807","20180824","20180913","20181001"];
subplotrow=4;
subplotcolumn=3;

for channel=6
    
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
    for sessionname=sessiondate
        % change directory to each day
        cd(strcat("/Volumes/Hippocampus/Data/picasso/",sessionname,"/",sessionno,""))
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
    
    for sessionname=sessiondate
        % Initialises subplot
        subplot(subplotrow,subplotcolumn,subplotno);
        
        % Generates subplot for spectrogram
        surf(chn.(strcat("Date_",sessionname,"")).T,chn.(strcat("Date_",sessionname,"")).F,chn.(strcat("Date_",sessionname,"")).Pnorm,"EdgeColor","none");
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
    saveas(gcf,strcat("/Users/",userfol,"/Dropbox/MATLAB/Hippocampal LFP Research/Theta/",channelname,"_spec.png",""))
    close
    
    % -------------------------------------------------------------------------
    % plot the non-standardised cc plot for all sessions with the same channel
    
    figure("Position", get(0, "Screensize"))
    subplotno=1;
    
    for sessionname=sessiondate
        % Initialises subplot
        subplot(subplotrow,subplotcolumn,subplotno);
        
        % Generates subplot for cc plot
        heatmap(chn.(strcat("Date_",sessionname,"")).F,flipud(chn.(strcat("Date_",sessionname,"")).F),flipud(chn.(strcat("Date_",sessionname,"")).cc),...
            "XLimits",{0,148.43750},"YLimits",{148.43750,0},...
            "XDisplayLabels",(round(chn.(strcat("Date_",sessionname,"")).F)),"YDisplayLabels",(flipud(round(chn.(strcat("Date_",sessionname,"")).F))),...
            "Colormap", jet,"ColorLimits",[-0.2,0.7],"ColorbarVisible","on",...
            "Title",strcat("session",string(sessionname),""),"FontSize",5);
        subplotno=subplotno+1;
    end
    
    sgtitle(strcat("CC plot averaged for all trials in each session for ",channelname,""),"FontSize",12);
    
    % Save the figure
    saveas(gcf,strcat("/Users/",userfol,"/Dropbox/MATLAB/Hippocampal LFP Research/Theta/",channelname,"_cc.png",""))
    close
    
    % Delete the structure array
    clear chn;
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

