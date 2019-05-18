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
    
    %     Spectrogram data for inter-trial interval
    idx = (tIdx(1)-(1000/1000*1000)):tIdx(1); %Inter-trial interval data (1000ms)
    data = vl.data.analogData(idx);
    datam = mean(data);
    [S,F,T,P]=spectrogram(data-datam,200,150,256,1000,'yaxis');
    
    %     Normalization parameters
    Pmean=mean(P,2); %mean power density of each frequency bin
    Pstd=std(P,0,2); %standard deviation of each frequency bin
    Smean=mean(S,2); %mean power density of each frequency bin
    Sstd=std(S,0,2); %standard deviation of each frequency bin
    
    %     Spectrogram data for trials
    idx = tIdx(1):tIdx(3); %Trial data
    data = vl.data.analogData(idx);
    datam = mean(data);
    [spec(n).S,spec(n).F,spec(n).T,spec(n).P,spec(n).Fc,spec(n).Tc]=...
        spectrogram(data-datam,200,150,256,1000,'yaxis');
    
    % trial psd and spectrogram normalised to the ITI period
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

% Normalize the time variable for the spectrogram information of the
% rest of the trial 

% Initialise the time-normalised matrix
normP=zeros(size(F,1),199,vl.data.numSets);
normS=zeros(size(F,1),199,vl.data.numSets);

for n=1:vl.data.numSets %for all trials
    %store the PSD and time data from 1s to end of trial to a temporary matrix
    tempNormT=spec(n).T(:,20:size(spec(n).T,2));
    tempNormP=spec(n).Pnorm(:,20:size(spec(n).Pnorm,2));
    tempNormS=spec(n).Snorm(:,20:size(spec(n).Snorm,2));
     
    %discretize the time stamps to fixed 0.050s widths (199 uniform
    %time bins)
    [binno,~] = discretize(tempNormT,199); %Create 199 uniform time bins
    tempno=1; %initialize normalized time marker
    
    for timeno=1:199
        if tempno==size(tempNormT,2) %If this is the last tempno, terminate
            normP(:,timeno,n)=tempNormP(:,tempno);
            normS(:,timeno,n)=tempNormS(:,tempno);
            continue
        elseif binno(tempno)==timeno 
            %If the bin matches the time number, add 1 to tempno
            if binno(tempno+1)==timeno+1 %If this bin is unique
                normP(:,timeno,n)=tempNormP(:,tempno);
                normS(:,timeno,n)=tempNormS(:,tempno);
            elseif binno(tempno+1)==timeno %If this and the next bin is identical
                normP(:,timeno,n)=(tempNormP(:,tempno)+tempNormP(:,tempno+1))/2; %Store the mean of the two powers
                normS(:,timeno,n)=(tempNormS(:,tempno)+tempNormS(:,tempno+1))/2;
            else %If this and the next bin is not identical or 1 different (next time bin has no power as it is stretched out)
                normP(:,timeno,n)=tempNormP(:,tempno); %this and the next time bin will have the same power
                normS(:,timeno,n)=tempNormS(:,tempno);
            end
            tempno=tempno+1;
        else %If the bin does not match the time number, do not add 1 to tempno
            normP(:,timeno,n)=tempNormP(:,tempno); %this timeno will have the same power as the last tempno
            normS(:,timeno,n)=tempNormS(:,tempno);
        end
    end
end 

% =========================================================================

% Creates a structure that averages the data across all trials

mcuespecP=zeros(size(F,1),218,vl.data.numSets); % Initialise variables
mcuespecS=zeros(size(F,1),218,vl.data.numSets); % Initialise variables
mspec=struct(); % Initialise variable

for n=1:vl.data.numSets
    % append the PSD information of the first 1s during cue
    % presentation to one big matrix
    mcuespecP(:,1:19,n)=spec(n).Pnorm(:,1:19);
    % append the normalised PSD information of the rest of the trial to one big matrix
    mcuespecP(:,20:218,n)=normP(:,:,n);
    
    % repeat for spectrogram information
    mcuespecS(:,1:19,n)=spec(n).Snorm(:,1:19);
    mcuespecS(:,20:218,n)=normS(:,:,n);
end

% Average trials that have the same cue during cue presentation
mspec.Pnorm=mean(mcuespecP,3);
mspec.Snorm=mean(mcuespecS,3);
% Recreate the frequency and time marker data
mspec.F=spec.F;
mspec.T=[0.1:0.05:10.95];

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

% =========================================================================

% Creates a structure that averages the normalised data across all trials
% without standardising the time bins

% Initialise variables
sizePnorm=zeros(vl.data.numSets);

% Obtain the largest size of Pnorm time data
for n=1:vl.data.numSets
    sizePnorm(n)=size(spec(n).Pnorm,2);
end

for s=1:max(sizePnorm)
    marker=1;
    Pbigmatrix=[];
    Sbigmatrix=[];
    for n=1:vl.data.numSets
        if s>size(spec(n).Pnorm,2)
            continue
        elseif s==size(spec(n).Pnorm,2)
            mspec.T2=spec(n).T;
            Pbigmatrix(:,marker)=spec(n).Pnorm(:,s);
            Sbigmatrix(:,marker)=spec(n).Snorm(:,s);
            marker=marker+1;
        else
            % append the PSD and spectrogram information to one big matrix
            Pbigmatrix(:,marker)=spec(n).Pnorm(:,s);
            Sbigmatrix(:,marker)=spec(n).Snorm(:,s);
            marker=marker+1;
        end
    end
    % Average the data from all trials together
    mspec.Pnorm2(:,s)=mean(Pbigmatrix,2);
    mspec.Snorm2(:,s)=mean(Sbigmatrix,2);
end

% Compute the correlation coefficients between frequencies for each
% averaged spectrogram data
Pmean=mean(mspec.Pnorm2,2); %mean power density of each frequency spectrum
Pdiff=zeros((size(mspec.Pnorm2,1)),(size(mspec.Pnorm2,2))); % Initialises Pdiff

for r=1:size(mspec.Pnorm2,1)
    % Difference at each time point with the mean power density for that frequency
    Pdiff(r,:)=mspec.Pnorm2(r,:)-Pmean(r);
end

for r=1:size(mspec.Pnorm2,1)
    for c=1:size(mspec.Pnorm2,1)
        %Correlation coefficient for each frequency band
        mspec.cc2(r,c)=sum(Pdiff(r,:).*Pdiff(c,:))/sqrt(sum(Pdiff(r,:).*Pdiff(r,:),2)*sum(Pdiff(c,:).*Pdiff(c,:),2));
    end
end

