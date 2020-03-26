
% LFP script eyeses to run the functions rpl parallel and rpllfp,
% spectrogram plot, and power spectrum across different channels for a
% specific sessioneye dataset

% The script will also compute the vmlfp data normalised to the sessioneyes
% data

userfol="/Users/AlizarinMoon";
savefolder="/seseye/";

for date="20181001"
    session="sessioneye";
    
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
                    [mrawspecdat(channelno).data,mspecdat(channelno).data,mvlspecdat(channelno).data]=combfun(channel,array,session,date,userfol,savefolder);
                end
                
            case 2 % array 2
                for channelno=33:64
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    [mrawspecdat(channelno).data,mspecdat(channelno).data,mvlspecdat(channelno).data]=combfun(channel,array,session,date,userfol,savefolder);
                end
                
            case 3 % array 3
                for channelno=65:96
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    [mrawspecdat(channelno).data,mspecdat(channelno).data,mvlspecdat(channelno).data]=combfun(channel,array,session,date,userfol,savefolder);
                end
                
            case 4 % array 4
                for channelno=97:128
                    channel=channelname(channelno);
                    %saves the individual plots for each channel
                    [mrawspecdat(channelno).data,mspecdat(channelno).data,mvlspecdat(channelno).data]=combfun(channel,array,session,date,userfol,savefolder);
                end
        end
    end
end

save(strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,"/",date,"_seseye_data.mat"),"mrawspecdat","mspecdat","mvlspecdat")
    
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

function [mrawspec,mspec,mvlspec] = combfun(channel,array,session,date,userfol,savefolder)

% Initialises directory name
folder=strcat("/Volumes/Hippocampus/Data/picasso/",date,"/",session,"/",array,"/",channel);

if exist(folder, "dir")==7 % Check that the folder exists
    
    % define parameters
    c1=-3; c2=3;
    ps1=-0.5; ps2=2 ;
    
    % creates a folder for the channel
    cd(strcat(userfol,'/Dropbox/MATLAB'))
    mkdir(strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel))
    
    % Generate the lfp time-domain data (Frequency 0-150Hz)
    cd(folder)
    rp = rplparallel('auto');
    lfp = rpllfp('auto');
    
    % ===============================================================
    %plot the signal average
%     signalaverage(lfp, [2 10], [-25 25])
%     saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel,'/signalaverage_',date,'_',channel,'.png'))
%     close
%     
    % ===============================================================
    % plot the raw spectrogram
    
    % Voltage data should be retrieved for each trial with 200ms before the cue time (baseline data)
    T=1;
    
    for n = 1:length(rp.data.markers) % trial marker
        
        if rp.data.markers(n) == 6
        
        % Obtain trial index for the marker
        % Each datapoint increases by 0.0010s, starting from index 0-
        % so each time stamp corresponds to index using the following formula:
        % index = time*1000+1
        
        % Trial for eye data will be from marker 6 to marker 32
        % Alternative: The 200ms prior to the marker 6 (cue)
        tIdx(1) = round(rp.data.timeStamps(n)*1000+1);
        tIdx(2) = round(rp.data.timeStamps(n+2)*1000+1);
        
        % Spectrogram data for the entire trial including a 1s pre-trial period
        idx = (tIdx(1)-200:(tIdx(1)));
        data = lfp.data.analogData(idx);
        datam = mean(data);
        [~,spec(T).F,spec(T).T,spec(T).P]=spectrogram(data-datam,40,30,(0:50),1000,'yaxis');
        
        % Normalisation to Baseline Period:
        % This is the old methodology where the spectrogram data for each trial is normalised to a pre-specified baseline period within the trial
        % each trial is normalised to a baseline period then averaged
        
        % Spectrogram data for the normalisation period
        idx = (tIdx(1)-500:tIdx(1)-200);
        data = lfp.data.analogData(idx);
        datam = mean(data);
        [~,~,~,P]=spectrogram(data-datam,200,150,(0:50),1000,'yaxis');
        
        % Normalisation parameters using the mean and stdev of the power
        % for the normalisation period
        Pmean=mean(P,2); %mean power density of each frequency bin
        Pstd=std(P,0,2); %standard deviation of each frequency bin
        
        spec(T).Pnorm=(spec(T).P-Pmean)./Pstd;
        
        T=T+1;
        
        end
        
    end
    
    % ===============================================================
    % Creates a structure that averages the normalised data across all trials
    % without standardising the time bins (but limit the data to 2s only)
    
    for s=1:17 % limit the data from 0s to 2s only (average mean length of the trials)
        % initialise variables
        Pmatrix=[];
        for n=1:length(spec) % for each trial
            if s>size(spec(n).P,2) % skips the trials that have already ended
                continue
            else
                % append the PSD and spectrogram information to one big matrix
                Pmatrix(:,n)=spec(n).P(:,s);
                Pmatrixnorm(:,n)=spec(n).Pnorm(:,s);
            end
        end
        % Average the data from all trials together
        mrawspec.Pnorm(:,s)=mean(Pmatrix,2);
        mspec.Pnorm(:,s)=mean(Pmatrixnorm,2);
    end
    
    %Initialise the time and frequency data
    mrawspec.T=spec.T;
    mrawspec.F=spec.F;
    mspec.T=spec.T;
    mspec.F=spec.F;

    % ===============================================================
    
    plotspec(mrawspec,0.02,0.18,1,50,-5,5);
    % save the raw spectrogram plot
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel,'/PSD_rawspec_',date,'_',channel,'.png'))
    close
    
    % compute and save normalised spectrogram plot
    plotspec(mspec,0.02,0.18,1,50,-5,5);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel,'/norm_-500to-200_PSD_',date,'_',channel,'.png'))
    close
    
    % comput power spectrum
    Pnorm=mean(mspec.Pnorm,2);
    plot(mspec.F,Pnorm,'color','black','LineWidth',2)
    axis([1 50 -5 5]);
    set(gca,'FontSize',15); xticks(0:5:50); yticks(-5:2.5:5);
    xlabel('Frequency (Hx)'); ylabel('Z-scored Power');
    
    % plot reference lines
    hold on
    line([2 2],[-100 100])
    line([5 5],[-100 100])
    line([6 6],[-100 100])
    line([9 9],[-100 100])
    line([0 30],[0 0],'LineStyle','--')
    hold off
    
    % save power spectrum plot
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel,'/norm_-500to-200_PS_',date,'_',channel,'.png'))
    close
    
    
    % ===============================================================
    % ===============================================================
    % ===============================================================
    % Computing vmlfp data normalised to the seseyes data
    
    cd(strcat("/Volumes/Hippocampus/Data/picasso/",date,"/session01/",array,"/",channel))
    vl = vmlfp("auto");
    
    % Normalisation parameters using the mean and stdev of the power
    % for the normalisation period
    Pmean=mean(mrawspec.Pnorm,2); %mean power density of each frequency bin
    Pstd=std(mrawspec.Pnorm,0,2); %standard deviation of each frequency bin
    
    % Spectrogram data for the entire trial including a 1s pre-trial period
    for n = 1:vl.data.numSets % trial number
        
        % Spectrogram data for the entire trial including a 1s pre-trial period
        tIdx = vl.data.trialIndices(n,:); % Obtain trial indexes
        idx = (tIdx(1)-(1000/1000*1000)):tIdx(3);
        data = vl.data.analogData(idx);
        datam = mean(data);
        [~,vlspec(n).F,vlspec(n).T,vlspec(n).P]=spectrogram(data-datam,200,150,(0:50),1000,'yaxis');
        
        % compute normalised spectrogram
        vlspec(n).Pnorm=(vlspec(n).P-Pmean)./Pstd;
        
    end
    
    % ===============================================================
    % Creates a structure that averages the normalised data across all trials
    % without standardising the time bins (but limit the data to 7s only)
    
    for s=1:159 % limit the data from -1.0s to 7s only (average mean length of the trials)
        Pmatrix=[];
        for n=1:vl.data.numSets % for each trial
            if s>size(vlspec(n).Pnorm,2) % skips the trials that have already ended
                continue
            else
                % append the PSD and spectrogram information to one big matrix
                Pmatrix(:,n)=vlspec(n).Pnorm(:,s);
            end
        end
        % Average the data from all trials together
        mvlspec.Pnorm(:,s)=mean(Pmatrix,2);
    end
    
    %Initialise the time and frequency data
    mvlspec.T=(-0.9:0.05:7);
    mvlspec.F=spec.F;
    
    % save normalised spectrogram plot
    plotspec(mvlspec,-0.9,7,1,50,-3,3);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel,'/PSD_vlnormspec_',date,'_',channel,'.png'))
    close
    
    %power spectrum
    powerspectrum(mvlspec,1,50,-5,5);
    saveas(gcf,strcat(userfol,'/Dropbox/MATLAB/Hippocampal LFP Research/',date,savefolder,channel,'/norm_seseye_PS_',date,'_',channel,'.png'))
    close
    
else
    mrawspec.Pnorm=[]; mrawspec.F=[]; mrawspec.T=[];
    mspec.Pnorm=[]; mspec.F=[]; mspec.T=[];
    mvlspec.Pnorm=[]; mvlspec.F=[]; mvlspec.T=[];
end

end


