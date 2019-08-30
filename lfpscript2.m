
% LFP script 2 to run the functions vmlfp, signal average, spectrogram
% plot, and power spectrum across different directories/ sessions

userfol="/Volumes/Users/mervyn";

for date="20181001"
    session="session01";
    
    for arrayno=1:4
        % Initialises array name
        array=strcat("array0",string(arrayno),"");
        
        % Runs the conglomerated function for each channel in each array;
        % switch for speed
        switch arrayno
            case 1 % array 1
                for channelno=1:32
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    mspecdat(channelno)=combfun(channel,array,session,date,userfol);
                end
                
            case 2 % array 2
                for channelno=33:64
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    mspecdat(channelno)=combfun(channel,array,session,date,userfol);
                end
                
            case 3 % array 3
                for channelno=65:96
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    mspecdat(channelno)=combfun(channel,array,session,date,userfol);
                end
                
            case 4 % array 4
                for channelno=97:128
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    mspecdat(channelno)=combfun(channel,array,session,date,userfol);
                end
        end
    end
end

save(strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,"/",date,"_mspecdat.mat"),"mspecdat")
    
% -------------------------------------------------------------------------

function [channel] = channelname(channelno)

% Initialises channel name
if numel(num2str(channelno))==1
    channel=strcat("channel00",string(channelno),"");
elseif numel(num2str(channelno))==2
    channel=strcat("channel0",string(channelno),"");
elseif numel(num2str(channelno))==3
    channel=strcat("channel",string(channelno),"");
end

end

function [mspec1] = combfun(channel,array,session,date,userfol)

% Initialises directory name
folder=strcat("/Volumes/Hippocampus/Data/picasso/",date,"/",session,"/",array,"/",channel);

if exist(folder, "dir")==7 % Check that the folder exists
    
    % define parameters
    c1=-3; c2=3;
    ps1=-0.5; ps2=2 ;
    
    % creates a folder for the channel
    cd(strcat(userfol,'/Dropbox/MATLAB'))
    mkdir(strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel))
    
    % Generate the lfp time-domain data (Frequency 0-150Hz)
    cd(folder)
    vl=vmlfp('auto');
    
    % ===============================================================
    %plot the signal average
    signalaverage(vl, [2 10], [-25 25])
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/signalaverage_',date,'_',channel,'.png'))
    close
    
    % ===============================================================
    %plot the raw spectrogram
    
    % Voltage data should be retrieved for each trial with 1000ms before the cue time (baseline data)
    for n = 1:vl.data.numSets % trial number
        
        tIdx = vl.data.trialIndices(n,:); % Obtain trial indexes
        
        % Spectrogram data for the entire trial including a 1s pre-trial period
        idx = (tIdx(1)-(1000/1000*1000)):tIdx(3);
        data = vl.data.analogData(idx);
        datam = mean(data);
        [~,spec(n).F,spec(n).T,spec(n).P]=spectrogram(data-datam,200,150,(0:50),1000,'yaxis');
    end
    
    % ===============================================================
    % Creates a structure that averages the normalised data across all trials
    % without standardising the time bins (but limit the data to 7s only)
    
    for s=1:159 % limit the data from -1.0s to 7s only (average mean length of the trials)
        marker=1; % initialise variables
        Pmatrix=[];
        for n=1:vl.data.numSets % for each trial
            if s>size(spec(n).P,2) % skips the trials that have already ended
                continue
            else
                % append the PSD and spectrogram information to one big matrix
                Pmatrix(:,marker)=log(spec(n).P(:,s));
                marker=marker+1;
            end
        end
        % Average the data from all trials together
        mrawspec.Pnorm(:,s)=mean(Pmatrix,2);
    end
    
    %Initialise the time and frequency data
    mrawspec.T=(-0.9:0.05:7);
    mrawspec.F=spec.F;
    
    % ===============================================================
    
    plotspec(mrawspec,-0.9,7,1,50,-3,3);
    % save the raw spectrogram plot
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/PSD_rawspec_',date,'_',channel,'.png'))
    close
    
    % ===============================================================
    % plot the baseline normalised spectrogram
    
    % each trial is normalised to a baseline period then averaged
    mspec1 = mnormspec(vl,"normal","norm",[-1000 -501]);
    %spectrogram plot
    plotspec(mspec1,-0.9,7,1,50,c1,c2);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/BLnorm_-1000to-501_PSD_',date,'_',channel,'.png'))
    close
    %power spectrum
    [mspec1.Pnorm1, mspec1.Pnorm2, mspec1.Pnorm3, mspec1.Pnorm4]=powerspectrum(mspec1,1,30,ps1,ps2);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/BLnorm_-1000to-501_PS_',date,'_',channel,'.png'))
    close
    
    % ===============================================================
    % plot the baseline corrected, log normalised spectrogram
    mspec=mnormspec(vl,"log");
    plotspec(mspec,-0.9,7,1,50,c1,c2);
    % save the raw spectrogram plot
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/PSD_normspec',date,'_',channel,'.png'))
    close
    
    % baseline correction from period -500 to -1
    mspec=mnormspec(vl,"log","baseline",[-500 -1]);
    % Plot the spectrogram
    plotspec(mspec,-0.9,7,1,50,c1,c2);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/baseline_-500to-1_PSD_',date,'_',channel,'.png'))
    close
    % Average the power spectrum across defined epochs
    powerspectrum(mspec,1,30,ps1,ps2);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,'/',channel,'/baseline_-500to-1_PS_',date,'_',channel,'.png'))
    close
    
else
    mspec1.Pnorm=[]; mspec1.F=[]; mspec1.T=[];
    mspec1.Pnorm1=[]; mspec1.Pnorm2=[]; mspec1.Pnorm3=[]; mspec1.Pnorm4=[];
end

end


