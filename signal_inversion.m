% LFP script to view the signal average plot for one channel across
% different sessions at the same time to look for signal inversion

% Dropbox user folder and file name to save
userfol="/Users/AlizarinMoon";
savename="_SAnav_em.png";

% Initialises session, channels, subplot variables
sessionno="session01";
arrayname="array01";
sessiondate=["20180221","20180312","20180323","20180402","20180404","20180416","20180625","20180704","20180807","20180824","20180913","20181001"];
subplotrow=4;
subplotcolumn=3;
subplotindex=[1,4,7,10,2,5,8,11,3,6,9,12];
epoch=[0 1000];
ylimits=[-30 30];
bandwidth=[2 10];
subplotno=1;

for channel=6
    
    % Initialises channel name
    if numel(num2str(channel))==1
        channelname=strcat("channel00",string(channel),"");
    elseif numel(num2str(channel))==2
        channelname=strcat("channel0",string(channel),"");
    elseif numel(num2str(channel))==3
        channelname=strcat("channel",string(channel),"");
    end
    
    figure("Position", get(0, "Screensize"));

    % Creates a matrix which stores the relevant data for the channel in
    % each different session
    for sessionname=sessiondate
        % change directory to each day
        cd(strcat("/Volumes/Hippocampus/Data/picasso/",sessionname,"/",sessionno,"/",arrayname,"/",channelname,""))
        vl = vmlfp('auto'); %computes the LFP data
        
        for n=1:vl.data.numSets
            tIdx = vl.data.trialIndices(n,:); %Obtain trial indexes
            idx = (tIdx(2)+epoch(1)):(tIdx(2)+epoch(2)); %Retrieve the epoch numbers (sampling rate is 1000)
            data = vl.data.analogData(idx); %retrieve data for the trial
            datam = mean(data);
            datamat(n,:) = data-datam; %save the demeaned data in a matrix
        end
        
        subplot(subplotrow,subplotcolumn,subplotno);
        dataaverage=mean(datamat);
        plot((epoch(1):epoch(2)),bandpass(dataaverage,bandwidth,1000),'.-') %signal average
        ylim(ylimits)
        title(sessionname,'FontSize',16)
        
        hold on
        line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
        hold off
        
        subplotno=subplotno+1;
    end
    
    sgtitle(strcat(channelname,"; Frequency range=",mat2str(bandwidth),"; Epoch referenced to start of navigation/ ms"),'FontSize',20)
    
    % Save the figure
    saveas(gcf,strcat(userfol,"/Dropbox/MATLAB/Hippocampal LFP Research/Theta/",channelname,savename,""))
    close
end
