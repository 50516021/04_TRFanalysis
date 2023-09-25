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
%%% 20230111 for Biosemi with 16 chnannel
%%% 20230318 adapted subject numbers
%%% v2 
%%% 20230322 adapted all BDF/EDF/CSV
%%% 20230328 proccessed data save
%%% 20230804 appried BPF [1.5 8], added date indentification

%%%%%%%%%!!check!!%%%%%%%%%%%
%%% name, tgArray, fpass

clearvars; close all;
% pname = 'R:\06_Tutorials\EEGanalysisWorkshop2022\example_ISNT';

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Chose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 

%%% get filenames
fname= strcat(outfolder, 'step3_epochs', '_', experiment_name, '_ICAprocessedAfterRejections.mat'); %Masker based epoched EEG data file name with its path
titlename = strrep(experiment_name, '_', ' '); %name for figure title (replased '_' to blank)
disp(['Processing: ' experiment_name]) %make sure the processing data

%%% get meta data (stimuli info)
metadata_file = ls([outfolder '/metadata*']);
load(metadata_file(1:end-1)); %participant's responces

% filter settings
fpass = [1.5 8]; %frequency of low/hi cut (Hz)
fsFilt = 230; %order of filtering

% trial information
Cols = ["blue", "red", "white", "green"]; %color variations
Nums = 1:8; %number variations
SNRs = [-12 -18]; %signal noise ratio variations
Spat = ["front", "back"]; %spatial pattern variations

numCols = length(Cols); %number of color variety
numNums = length(Nums); %number of number variety
numSNRs = length(SNRs); %number of SNR variety
numSpat = length(Spat); %number of Spatial pattern variety

baselinedur = 0.3; %duration of baseline (sec)
stimulidur = 32; %minimum target start time (sec)
targetdur   = 2.8; %target duration (sec)

%% trigger values

trg_grdst  = 250; %grand start/end
trg_instr  = 240; %instruction onset
trg_onset  = 230; %stream onset
trg_offset = 220; %stream offset

%% load file
load(fname) %load Masker based EEG
numChMsk = size(epochs,2); %number of channels on masker
epochs_temp=epochs; %temporal signal for filtering
for i =1:size(epochs,3)
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

t = -baselinedur: 1/fs : stimulidur; %sample to second conversion
streamdur = 256*(baselinedur+stimulidur); %plot stream duration

numTrial = size(epochs,3);

plotEp = epochs(:,[Hotch(1) Hotch(2)],:)-epochs(:,[Coldch(1) Coldch(2)],:)/2-epochs(:,[Coldch(2) Coldch(1)],:)/2; %subtraction from center channel
%index: plotEp(EEG, channel, trial)
saveEp = plotEp(round(baselinedur*fs)+1:end,:,:);

% for i = 1:numTrial
%     subplot(numTrial,1,i);
%     plot(t,squeeze(plotEp(:,:,i))); hold on;  
%     ylim([-yscale yscale]); 
%     xlim([-baselinedur stimulidur]); 
%     line([0;0], get(gca, 'ylim')); hold off;
%     % title(sprintf('Spatial Pattern: %s, %d trials',Spat(i), numTrial));
%     xlabel('time[s]');    
% end
% 
% sgtitle(sprintf('Final evoked responce [%s] for Masker onset BPF[%0.1f-%0.1f]', titlename, fpass))
% legend(legends(:), 'location', 'southeast');
% saveas(gcf, strcat(outfolder, 'final_EvokedResponce_Msk_', experiment_name, '_filt', '.pdf'))

%% save proccessed data
fs_EEG = fs;
date = datestr(now,'yyyymmdd');
save_filename = strcat(outfolder, 'step4_plotdata_', date, '_',  experiment_name, '.mat');
save(save_filename,'saveEp','metadata','fs_EEG');
disp(strcat(save_filename, ' has been saved'))


