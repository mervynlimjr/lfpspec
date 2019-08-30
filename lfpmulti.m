% LFP script to view multiple objects for individual trials of a channel in
% one session at the same time by loading instantiated and initialising relevant
% objects

% Initialises session, channels, subplot variables
sessionno="session01";
arrayname="array02";
sessiondate={'20181017','20180917','20180807','20180625'};
channel=43;

% Initialises channel name
if numel(num2str(channel))==1
    channelname=strcat("channel00",string(channel),"");
elseif numel(num2str(channel))==2
    channelname=strcat("channel0",string(channel),"");
elseif numel(num2str(channel))==3
    channelname=strcat("channel",string(channel),"");
end

% -------------------------------------------------------------------------
% Creates a new structure which stores the relevant data for the channel in
% each different session
% 
% for sessionname=sessiondate
%     % change directory to each day
%     cd (strcat('/Volumes/Hippocampus/Data/picasso/',sessionname,'/',sessionno,""));
%     % load pre-instantiated um
%     load("unitymaze.mat");
%     % retrieve the data for the required channel
%     chn.(strcat("Date_",sessionname,"")).um=um;
%     % clear existing data
%     clear um;
%     
%     % change directory to the channel directory
%     cd (strcat('/Volumes/Hippocampus/Data/picasso/',sessionname,'/',sessionno,'/',arrayname,"/",channelname,""));
%     % generate vr, vh, and vl
%     chn.(strcat("Date_",sessionname,"")).vl=vmlfp('auto');
%     chn.(strcat("Date_",sessionname,"")).vr=vmhighpass('auto','Raw');
%     chn.(strcat("Date_",sessionname,"")).vh=vmhighpass('auto');
% end

% -------------------------------------------------------------------------
% Brute force through InspectGUI to look at up to 4 different channels/
% days

InspectGUI(chn.(strcat("Date_",sessiondate(1),"")).um,'addObjs',...
    {chn.(strcat("Date_",sessiondate(2),"")).um,chn.(strcat("Date_",sessiondate(3),"")).um,chn.(strcat("Date_",sessiondate(4),"")).um,...
    chn.(strcat("Date_",sessiondate(1),"")).vr,chn.(strcat("Date_",sessiondate(2),"")).vr,chn.(strcat("Date_",sessiondate(3),"")).vr,chn.(strcat("Date_",sessiondate(4),"")).vr,...
    chn.(strcat("Date_",sessiondate(1),"")).vh,chn.(strcat("Date_",sessiondate(2),"")).vh,chn.(strcat("Date_",sessiondate(3),"")).vh,chn.(strcat("Date_",sessiondate(4),"")).vh,...
    chn.(strcat("Date_",sessiondate(1),"")).vl,chn.(strcat("Date_",sessiondate(2),"")).vl,chn.(strcat("Date_",sessiondate(3),"")).vl,chn.(strcat("Date_",sessiondate(4),"")).vl,...
    chn.(strcat("Date_",sessiondate(1),"")).vl,chn.(strcat("Date_",sessiondate(2),"")).vl,chn.(strcat("Date_",sessiondate(3),"")).vl,chn.(strcat("Date_",sessiondate(4),"")).vl,...
    chn.(strcat("Date_",sessiondate(1),"")).vl,chn.(strcat("Date_",sessiondate(2),"")).vl,chn.(strcat("Date_",sessiondate(3),"")).vl,chn.(strcat("Date_",sessiondate(4),"")).vl},...
    'SP',[6 4],'optArgs',{...
    {'Trial'},{'Trial'},{'Trial'},{'Trial'},...
    {},{},{},{},...
    {},{},{},{},...
    {},{},{},{},...
    {'Filter','FilterWindow',[5 15]},{'Filter','FilterWindow',[5 15]},{'Filter','FilterWindow',[5 15]},{'Filter','FilterWindow',[5 15]}...
    {'TFfft','TFfftPoints',512},{'TFfft','TFfftPoints',512},{'TFfft','TFfftPoints',512},{'TFfft','TFfftPoints',512}});
