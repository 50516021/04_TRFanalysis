%%% mTRF EEG analysis for EDF/CSV - step 1 ver4 %%% 
%%% - epoching 
%%% - feasibility test 2023
%%%
%%% required Add-ons
%%% - EEGLAB
%%%
%%% required functions
%%% - BPF.m
%%% - openbdf.m
%%% - readbdf.m

%%% v1  
%%% 20230920 experiment 'experiment_mTRF_feasibility_v1.m'
%%% v2
%%% 20230923 experiment 'experiment_mTRF_feasibility_v2.m'
%%%     NEXT -> step4
%%% 20231010 OS flex
%%% 20231115 also for v3
%%% 20231129 for v4

addpath('../');
addpath('../../02_EEGanalysis'); %add path of EEGanalysis
addpath('../../02_EEGanalysis/eeglab2023.0'); %add path of EEGLAB

OSflag = OSdetection_v1;

%% parameter

%%% triggers %%%
trg_grdst  = 250;
trg_onset  = 230;
trg_instr  = [210, 200]; %instruction trigger
trg_offset = 220;

%%% duration = starttime + streamdur + poststream (refer experiment diagram)
pretrigger = 0.5;  %time duration before trigger (sec)
poststream = 0.5;  %time duration after stimuli (sec)

baselinedur= 0.3;  %duration of baseline (sec)
triggerthre = 0.70; %threshould of the trigger (sec)
% streamdur = 32.0; %stream duration (sec)
fpass = [1.5 20]; %frequency of low/hi cut (Hz)
fs = 256; %resample rate

%% experiment info 

%%% determine the device
devOpt = ["Biosemi (BDF)" "DSI-24 (CSV/EDF)"]; % EEG device option
prompt = 'Choose EEG device'; % prompt message
[dev,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',devOpt); % option selection window
datadir = {'../../01_OriginalData/BDFdata/' '../../01_OriginalData/EDFdata/'};

%%% subject number %%%
prompt = {'Enter subjects number:'}; 
dlgtitle = 'Input';
dims = [1 35];
definput = {'00000'};
answerNum = inputdlg(prompt,dlgtitle,dims,definput);
sNum = answerNum{1};

%%% experiment name %%%
folders = struct2table(dir([datadir{dev}, 's', sNum, '*']));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[250 200],'ListString',folders.name); % option selection window
experiment_name = folders.name(foldInd,:); %subject (experiment) name

if iscell(experiment_name)
    experiment_name = experiment_name{:};
end

outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 
savepath = '';

if ~exist(outfolder, 'dir') %check folder existance
    mkdir(outfolder) %make directory
end

%% get behavioral data

filepath = strcat(datadir{dev}, experiment_name, '/');
metadata_file = ls([filepath, 'metadata*']); %find metadata (stimuli info) file

if OSflag(1) == "1" %Mac
    metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
    load(metadata_file); 
    save([outfolder extractAfter(metadata_file, filepath)], 'inst_flg', 'path_Tgt', 'timerange_Tgt', 'path_Msk', 'timerange_Msk'); %copy metadata file 

elseif OSflag(1) == "2" %Windows
    load([filepath metadata_file]); 
    save([outfolder metadata_file], 'inst_flg', 'path_Tgt', 'timerange_Tgt', 'path_Msk', 'timerange_Msk'); %copy metadata file 

end

%% load EEG data 
%%% Biosemi (BDF) %%%
if dev == 1
    %%% BDF file name %%%
    files = struct2table(dir([datadir{dev}, 's', sNum, '*/*bdf']));
    prompt = 'Chose BDF file name:';  % prompt message
    [foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[250 200],'ListString',files.name); % option selection window
    if ~iscell(files.folder)
        bdfname = strcat(files.folder(foldInd,:), '/', files.name(foldInd,:)); %subject (experiment) name
    else
        bdfname = strcat(files.folder{foldInd,:}, '/', files.name{foldInd,:}); %subject (experiment) name
    end

    h = openbdf(bdfname); %open BDF file
    d = readbdf(h,1:h.Head.NRec); %read BDF data
    tg = d.Record(end,:)'; %trigger 

%     numCh = 16*floor(size(d.Record,1)/16); %number of EEG channels (except trigger)
    numCh = 16; %since Biosemi only records 32 channels minimum, the number of the channles is set manually
    fsOgnl = h.Head.SampleRate(1); %Original EEG sampling rate
    eeg = BPF(d.Record(1:numCh,:)', fsOgnl, fsOgnl, fpass(1), fpass(2)); %band pass filter

    onsets = find((fs*1>diff(tg)).*(diff(tg)>0)) + 1; %index of all peaks in tg channel
    events = tg(onsets) - tg(onsets-2); %trigger value

%     %%% save original EEG data %saving BDF file is very big (because the sampling rate is too big)
%     save([outfolder, 'original_', experiment_name '.mat'],'h','d');

%%% DSI-24 (CSV/EDF) %%%
elseif dev == 2
    eeglab %read DSI-24 data
    fprintf('Load CSV data from EEGLAB window or load the .mat file manually (Subject No.%s)\n', convertCharsToStrings(sNum));
    keyboard; %wait
    
    EEGstruct = EEG; %EEG data acquired by EEGLAB
    chConv = [1:15 18:19 21:23]; %channel extraction (9-A1,16-X3,17-X2,18-A2,20-X1)
    numCh = length(chConv); %number of EEG channels (except trigger)
    fsOgnl = EEGstruct.srate; %[manual]

    eeg = EEGstruct.data; %get EEG data
    eeg = eeg(chConv,:); %reorder channels
    eeg = BPF(double(eeg)', fsOgnl, fsOgnl, fpass(1), fpass(2)); % BPF = band pass filter

    %%% trigger data
    Ev = EEGstruct.event; %[manual] 
    onsets = [Ev.latency]; 
    onsets = onsets(:)'; %Hwan: 5-end, SY: 3-end [manual]
    events = {Ev.type};
    events = str2double(events(:))'; %Hwan: 5-end, SY: 3-end [manual]

    %%% save original EEG data
    save([outfolder, 'original_', experiment_name '.mat'],'EEG');
end

%%% get onset's indices
stimOnsets = onsets(events==trg_onset); %index of simuli onset 
streamdur = timerange_Tgt(2) - timerange_Tgt(1);
numTrial = length(stimOnsets);

%% baseline
t = -baselinedur:1/fs:(streamdur); %full time base
basetime = -baselinedur: 1/fsOgnl : (streamdur - 1/fsOgnl); %base time of each stream 
baseline = find( (-baselinedur<basetime) .* (basetime<0) ); %EEG baseline guide

%% epoching 
clear epochs %clear cash
streamdur_epoch = fix(fsOgnl*(baselinedur+streamdur))+1; %data stream duration
epochs    = zeros(fix(fs*(baselinedur+streamdur))+1,numCh,numTrial); %epochs(EEG signal, channel, trial)

for k = 1:numTrial
    % if stimOnsets(k)+fix(fsOgnl*streamdur) < size(eeg,1)
        temp = eeg(stimOnsets(k)-fix(fsOgnl*baselinedur)+1 : stimOnsets(k)+fix(fsOgnl*streamdur),:); %EEG data region following masker onset triggers
    % else
    %     temp = eeg(stimOnsets(k)-fix(fsOgnl*baselinedur)+1 : stimOnsets(k)+fix(fsOgnl*(streamdur-40)),:);
    % end
    temp = temp - repmat(mean(temp(baseline,:)),size(temp,1),1); %subtract baseline
    epochs(:,:,inst_flg(k)) = resample(temp,fs,fsOgnl);     %subtracted epoch
end

epochs = epochs(fs*baselinedur+1:end,:,:);
save([outfolder 'step1_' experiment_name '.mat'],'epochs','t','fs');

