function [] = signalaverage(vl, bandwidth, ylimits)
% SIGNALAVERAGE computes the average of the time domain data and displays
% it as a figure
%
%   Inputs: lfpfig(vl,bandwidth,epoch)
%   (1) vl = vmlfp object - raw lfp data
%   (2) bandwidth = [x y], where x and y specifies the frequency limits of
%   the signal to be displayed
%   (3) ylimits = [a b], where a and b specifies the y-limits of the fig
%
%   Outputs:
%   A figure computing the odd and even signal averages

% =========================================================================
% plot the signal averages

epoch=[-2000 -1000];

for n=1:vl.data.numSets
    tIdx = vl.data.trialIndices(n,:); %Obtain trial indexes
    idx = (tIdx(2)+epoch(1)):(tIdx(2)+epoch(2)); %Retrieve the epoch numbers (sampling rate is 1000)
    data = vl.data.analogData(idx); %retrieve data for the trial
    datam = mean(data);
    datamat(n,:) = data-datam; %save the demeaned data in a matrix
end

dataodd=datamat(1:2:end,:);
oddaverage=mean(dataodd);
dataeven=datamat(2:2:end,:);
evenaverage=mean(dataeven);
totalaverage=mean(datamat);
subtractaverage=(oddaverage-evenaverage)/2;

figure("Position", get(0, "Screensize"));
subplot(3,3,1);
plot((epoch(1):epoch(2)),bandpass(oddaverage,bandwidth,1000),'.-') %signal average of the odd results
hold on
plot((epoch(1):epoch(2)),bandpass(evenaverage,bandwidth,1000),'.-') %signal average of the even results

% Plot lines to mark the zero crossing
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
ylim(ylimits)
ylabel('Superimposed Odd (Blue) and (Orange)','FontSize',13)
hold off

subplot(3,3,4);
plot((epoch(1):epoch(2)),bandpass(totalaverage,bandwidth,1000),'.-') %total average
ylim(ylimits)
ylabel('Signal Average = Odd+Even','FontSize',13)
% Plot lines to mark the zero crossing
hold on
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
hold off

subplot(3,3,7);
plot((epoch(1):epoch(2)),bandpass(subtractaverage,bandwidth,1000),'.-') %subtraction of odd and even
ylim(ylimits)
ylabel('+/: Signal Average = Odd-Even','FontSize',13)
% Plot lines to mark the zero crossing
hold on
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
hold off

% =========================================================================

epoch=[-1000 0];
for n=1:vl.data.numSets
    tIdx = vl.data.trialIndices(n,:); %Obtain trial indexes
    idx = (tIdx(2)+epoch(1)):(tIdx(2)+epoch(2)); %Retrieve the epoch numbers (sampling rate is 1000)
    data = vl.data.analogData(idx); %retrieve data for the trial
    datam = mean(data);
    datamat(n,:) = data-datam; %save the demeaned data in a matrix
end

dataodd=datamat(1:2:end,:);
oddaverage=mean(dataodd);
dataeven=datamat(2:2:end,:);
evenaverage=mean(dataeven);
totalaverage=mean(datamat);
subtractaverage=(oddaverage-evenaverage)/2;

subplot(3,3,2);
plot((epoch(1):epoch(2)),bandpass(oddaverage,bandwidth,1000),'.-') %signal average of the odd results
hold on
plot((epoch(1):epoch(2)),bandpass(evenaverage,bandwidth,1000),'.-') %signal average of the even results

% Plot lines to mark the zero crossing
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
ylim(ylimits)
hold off

subplot(3,3,5);
plot((epoch(1):epoch(2)),bandpass(totalaverage,bandwidth,1000),'.-') %total average
ylim(ylimits)
% Plot lines to mark the zero crossing
hold on
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
hold off

subplot(3,3,8);
plot((epoch(1):epoch(2)),bandpass(subtractaverage,bandwidth,1000),'.-') %subtraction of odd and even
ylim(ylimits)
% Plot lines to mark the zero crossing
hold on
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
hold off

xlabel('Epoch referenced to start of navigation/ ms','FontSize',16)

% =========================================================================

epoch=[0 1000];
for n=1:vl.data.numSets
    tIdx = vl.data.trialIndices(n,:); %Obtain trial indexes
    idx = (tIdx(2)+epoch(1)):(tIdx(2)+epoch(2)); %Retrieve the epoch numbers (sampling rate is 1000)
    data = vl.data.analogData(idx); %retrieve data for the trial
    datam = mean(data);
    datamat(n,:) = data-datam; %save the demeaned data in a matrix
end

dataodd=datamat(1:2:end,:);
oddaverage=mean(dataodd);
dataeven=datamat(2:2:end,:);
evenaverage=mean(dataeven);
totalaverage=mean(datamat);
subtractaverage=(oddaverage-evenaverage)/2;

subplot(3,3,3);
plot((epoch(1):epoch(2)),bandpass(oddaverage,bandwidth,1000),'.-') %signal average of the odd results
hold on
plot((epoch(1):epoch(2)),bandpass(evenaverage,bandwidth,1000),'.-') %signal average of the even results

% Plot lines to mark the zero crossing
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
ylim(ylimits)
hold off

subplot(3,3,6);
plot((epoch(1):epoch(2)),bandpass(totalaverage,bandwidth,1000),'.-') %total average
ylim(ylimits)
% Plot lines to mark the zero crossing
hold on
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
hold off

subplot(3,3,9);
plot((epoch(1):epoch(2)),bandpass(subtractaverage,bandwidth,1000),'.-') %subtraction of odd and even
ylim(ylimits)
% Plot lines to mark the zero crossing
hold on
line([epoch(1) epoch(2)],[0 0],[10 10; 10 10],'Color','k');
hold off

sgtitle(strcat(pwd," (frequency range=",mat2str(bandwidth),")"),'FontSize',20)

% =========================================================================

end


