function [spec, mspec] = lfpspec(vl)
% LFPSPEC Computes the normalised spectrogram and correlation coefficient
% data from the time domain data of the LFPs average across all trials in
% one session
%
%   Inputs:
%   vl = vmlfp object where vl.data comprises of the following structures:
%   analogData, analogInfo, analogTime, numSets, Args, markers, timeStamps,
%   and trialIndices
%
%   Outputs:
%   spec is a structure containing the normalised spectrogram and
%   correlation coefficient data for each trial
%   mspec is a structure containing the averaged normalised spectrogram and
%   correlation coefficient data across all trials

% =========================================================================

spec(vl.data.numSets)=struct(); % Initialise the normalised spectrogram matrix

% Voltage data should be retrieved for each trial with 1000ms before the
% cue time (baseline data)
for n = 1:vl.data.numSets % trial number
    
    tIdx = vl.data.trialIndices(n,:); % Obtain trial indexes
    
    %     Spectrogram data for the 'normalisation period NP'
    idx = (tIdx(1)-(1000/1000*1000)):(tIdx(1)-(501/1000*1000)); %Inter-trial interval data (-1000ms to -500ms)
    data = vl.data.analogData(idx);
    datam = mean(data);
    [S,F,T,P]=spectrogram(data-datam,200,150,512,1000,'yaxis');
    
    %     Normalization parameters of the NP
    Pmean=mean(P,2); %mean power density of each frequency bin
    Pstd=std(P,0,2); %standard deviation of each frequency bin
    Smean=mean(S,2); %mean power density of each frequency bin
    Sstd=std(S,0,2); %standard deviation of each frequency bin
    
    %     Spectrogram data for trials
    idx = (tIdx(1)-(500/1000*1000)):tIdx(3); %Trial data including the pre-trial data of -500ms
    data = vl.data.analogData(idx);
    datam = mean(data);
    [spec(n).S,spec(n).F,spec(n).T,spec(n).P,spec(n).Fc,spec(n).Tc]=...
        spectrogram(data-datam,200,150,512,1000,'yaxis');
    
    % trial psd and spectrogram normalised to NP
    spec(n).Pnorm=(spec(n).P-Pmean)./Pstd;
    spec(n).Snorm=(spec(n).S-Smean)./Sstd;
    
    %     Computes correlation coefficients for each trial
    CCmean=mean(spec(n).Pnorm,2); %mean power density of each frequency spectrum of trials
    Pdiff=zeros((size(spec(n).Pnorm,1)),(size(spec(n).Pnorm,2))); % Initialises Pdiff

    for r=1:size(spec(n).Pnorm,1)
        % Difference at each time point with the mean power density for that frequency
        Pdiff(r,:)=spec(n).Pnorm(r,:)-CCmean(r);
    end
    
    for r=1:size(spec(n).Pnorm,1)
        for c=1:size(spec(n).Pnorm,1)
            %Correlation coefficient for each frequency band
            spec(n).cc(r,c)=sum(Pdiff(r,:).*Pdiff(c,:))/...
                sqrt(sum(Pdiff(r,:).*Pdiff(r,:),2)*sum(Pdiff(c,:).*Pdiff(c,:),2));
        end
    end
end

% =========================================================================

% Creates a structure that averages the normalised data across all trials
% without standardising the time bins (but limit the data to 11s only)

for s=1:151 % limit the data from -05s to 7s only (average mean length of the trials)
    marker=1; % initialise variables
    Pmatrix=[];
    Smatrix=[];
    for n=1:vl.data.numSets % for each trial
        if s>size(spec(n).Pnorm,2) % skips the trials that have already ended
            continue
        else
            % append the PSD and spectrogram information to one big matrix
            Pmatrix(:,marker)=spec(n).Pnorm(:,s);
            Smatrix(:,marker)=spec(n).Snorm(:,s);
            marker=marker+1;
        end
    end
    % Average the data from all trials together
    mspec.Pnorm(:,s)=mean(Pmatrix,2);
    mspec.Snorm(:,s)=mean(Smatrix,2);
end

%Initialise the time and frequency data
mspec.T=(-0.5:0.05:7);
mspec.F=spec.F;

% Compute the correlation coefficients between frequencies for each
% averaged spectrogram data
Pmean=mean(mspec.Pnorm,2); %mean power density of each frequency spectrum
Pdiff=zeros((size(mspec.Pnorm,1)),(size(mspec.Pnorm,2))); % Initialises Pdiff

for r=1:size(mspec.Pnorm,1)
    % Difference at each time point with the mean power density for that frequency
    Pdiff(r,:)=mspec.Pnorm(r,:)-Pmean(r);
end

for r=1:size(mspec.Pnorm,1)
    for c=1:size(mspec.Pnorm,1)
        %Correlation coefficient for each frequency band
        mspec.cc(r,c)=sum(Pdiff(r,:).*Pdiff(c,:))/sqrt(sum(Pdiff(r,:).*Pdiff(r,:),2)*sum(Pdiff(c,:).*Pdiff(c,:),2));
    end
end

