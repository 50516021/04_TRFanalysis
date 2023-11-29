%%% Spatial Attention EEG analysis for BDF/EDF/CSV - step 4 %%% 
%%% - evoked responce
%%%
%%% required Add-ons
%%% - EEGLAB
%%% - Symbolic Math Toolbox
%%% required functions
%%% - 
%%% required setting files
%%% - 

%%% v1  
%%% 20230920 experiment 'experiment_mTRF_feasibility_v1.m'
%%% v2
%%% 20230923 experiment 'experiment_mTRF_feasibility_v2.m'
%%% 20231010 OS flex
%%% 20231129 also for v4

clearvars; close all;

OSflag = OSdetection_v1;

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Chose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 

if OSflag(1) == "1"
    %%% get filenames
    EEGfile = ls([outfolder, 'step1_*']); %find responce file
    EEGfile = EEGfile(1:end-1); %extract unnecessary charactar
    load(EEGfile); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
    load(metadata_file); %participant's responces

elseif OSflag(1) == "2"
    %%% get filenames
    EEGfile = ls([outfolder, 'step1_*']); %find responce file
    load([outfolder EEGfile]); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    load([outfolder metadata_file]); %participant's responces
end

stimulidur = timerange_Tgt(2) - timerange_Tgt(1);

% filter settings
fpass = [1.5 8]; %frequency of low/hi cut (Hz)
fsFilt = 230; %order of filtering

%% prepare EEG
numChMsk = size(epochs,2); %number of channels on masker
epochs_temp=epochs; %temporal signal for filtering
numTrial = size(epochs,3);

for i =1:numTrial
    epochs(:,:,i) = BPF(double(squeeze(epochs_temp(:,:,i))), fs, fsFilt, fpass(1), fpass(2)); % BPF = band pass filter
end

%% make figure

% figure;
    % 4...Fz
    % 8...Cz
    % 12..Pz (DSI-24>No Pz)
    % 14..O1
    % 16..O2 (DSI-24>15..O2)
if size(epochs,2) == 20 %DSI-24
    disp('Device: DSI-24')
    Hotch = [4 8]; % Fz and Cz
    Coldch = [14 15]; % O1 and O2
else %Biosemi
    disp('Device: Biosemi')
    Hotch = [4 8]; % Fz and Cz
    Coldch = [14 16]; % O1 and O2
end

ch1 = 'Fz';
ch2 = 'Cz';

yscale = 10; %scale of y axis
legends = {ch1, ch2, 'Onset'};

t = 1: 1/fs : stimulidur; %sample to second conversion
streamdur = fs*stimulidur; %plot stream duration

plotEp = epochs(:,[Hotch(1) Hotch(2)],:)-epochs(:,[Coldch(1) Coldch(2)],:)/2-epochs(:,[Coldch(2) Coldch(1)],:)/2; %subtraction from center channel
%index: plotEp(EEG, channel, trial)
saveEp = plotEp;

%% save proccessed data
fs_EEG = fs;
date = datestr(now,'yyyymmdd');
save_filename = strcat(outfolder, 'step4_plotdata_', date, '_',  experiment_name, '.mat');
save(save_filename,'saveEp', 'epochs', 'fs_EEG');
disp(strcat(save_filename, ' has been saved'))


