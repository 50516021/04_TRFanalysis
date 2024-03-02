%%% mTRF analysis - step 5 mTRF model generation for various stimuli length - nonsub -> use raw data of all the cannels %%% 
%%% - evoked responce, use all channels, two instructions version
%%%
%%% required Add-ons
%%% - 
%%% - 
%%% required functions
%%% - StimTrans_mTRF_v1.m
%%% - TRFestimation_v1.m
%%%
%%% required setting files
%%% - 

%%% v1  
%%% 20230921 experiment 'experiment_mTRF_feasibility_v1.m'
%%% v2
%%% 20230923 experiment 'experiment_mTRF_feasibility_v2.m'
%%% v3
%%% 20230927 for signal qualities   
%%% 20231129 nonsub - raw data
%%% v4
%%% 20231207 two instructions 'experiment_mTRF_feasibility_v4.m'
%%% 20240222 ICA option

clearvars; 
close all;

addpath('../'); %add path above
addpath('../../02_EEGanalysis'); %add path of EEGanalysis

OSflag = OSdetection_v1;

ICAopt = 1; %use ICA data (1) or not (0)
if ICAopt
    namekey = 'step4_plotdatav3_*';
    nameopt = "_ICA_";
else      
    namekey = 'step4_plotdata_*';   
    nameopt = "_";
end

Hotch = [4 8]; % Fz and Cz
Coldch = [14 15]; % O1 and O2

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 

if ICAopt
    outfolder_mTRFfig = strcat(outfolder, 'mTRF_fig/');
    outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl/');
else
    outfolder_mTRFfig = strcat(outfolder, 'mTRF_fig_nonsub/');
    outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl_nonsub/');
end

mkdir(outfolder_mTRFfig)
mkdir(outfolder_mTRFmdl)

if OSflag(1) == "1" %Mac
    %%% get filenames
    EEGfile = ls([outfolder, namekey]); %find responce file
    EEGfile = EEGfile(1:end-1); %extract unnecessary charactar
    load(EEGfile); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
    load(metadata_file); %participant's responces

elseif OSflag(1) == "2" %Windows
    %%% get filenames
    EEGfile = ls([outfolder, namekey]); %find responce file
    load([outfolder EEGfile]); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    load([outfolder metadata_file]); %participant's responces
end

if ICAopt; epochs = saveEp; end %use subtracted (ICAed) epoch
% epochs(:,[Hotch(1) Hotch(2)],:) = epochs(:,[Hotch(1) Hotch(2)],:)-epochs(:,[Coldch(1) Coldch(2)],:)/2-epochs(:,[Coldch(2) Coldch(1)],:)/2; %subtraction from center channel


%% parameters

stim_tag = ["Left", "Right", "Mixed"]; 
% stim_dur = [5*60 10*60 stimulidur]; %stimuli extraction duration [sec] 

%% stimulus preparation

stimfilename = "stimulus_mTRF/stim_mTRfpilot_v4.mat";
ExpVer = contains(experiment_name,'mTRFpilot_v4'); %experiment version check
ExpStim = exist(stimfilename, 'file');

%%% stimuli info
fs_Sound_down = 8000; %consider the frequency range of speech 
stimulidur = timerange_Tgt(2) - timerange_Tgt(1);

if ExpVer && ExpStim
    disp('### loading stimuli ###')
    stim_str = load(stimfilename);
    stim = stim_str.stim;
    disp('### stimuli have been loaded ###')
else
    [stimulus, fs_Sound] = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk);
    stimulus(:,3) = stimulus(:,1) + stimulus(:,2); %mixed stimulus
    
    %%% resample and extract stimuli
    stimulus_downsmpl = resample(stimulus,fs_Sound_down,fs_Sound);
    
    for i =1:size(stimulus,2)
        disp('### begin stimuli conversion ###')
        stim(:,:,i) = StimTrans_mTRF_v1(stimulus_downsmpl(:,i), fs_Sound_down);
        disp('### conversion done ###')
    end
    
    filename_stim = strcat(outfolder, 'stim_', experiment_name, '.mat');
    save(filename_stim,'stim');
end

%% data preparation

%%% chanel info
chs = ["Fz", "Cz"];
numch = [4, 8];
instruction = ["Left", "Right"];

fs_New = 128;

%% mTRf processing

for i = 1:length(inst_flg)
    for j =1:size(stim,3)
        stim_ext = resample(stim, fs_New, fs_Sound_down);
        resp = resample(epochs(1:stimulidur*fs_EEG,:,i), fs_New, fs_EEG);   
        
        %%% model training
        model = TRFestimation_v1(stim_ext(:,:,j), fs_New, resp, fs_New, 0, 0);
        filename = sprintf('mTRF_%s_inst%s', stim_tag(j), instruction(i));
        filename_mdl = strcat(outfolder_mTRFmdl, 'model', nameopt,  filename, '.mat');
        save(filename_mdl,'model');

        %%% mTRF estimation (figure)
        for k = 1:length(chs)     
            
            TRFestimation_v1(stim_ext(:,:,j), fs_New, resp, fs_New, k, model);

            sgtitle(sprintf('mTRF stimulus: %s,inst: %s, %s ', stim_tag(j), instruction(i), chs(k)))
            filename_pdf = strcat(outfolder_mTRFfig, filename, nameopt, chs(k), '.pdf');
            saveas(gcf, filename_pdf)
            disp(strcat(filename, ' figure has been saved'))

        end    
    end
end

%% stimulus preparation function %%%
function [stimulus, fs_t] = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk)
    %%% 2023213 resample section

    d_fifo = 10; %duration of fadein/fadeout (ms)
    
    disp('preparing stimuli')
    [wavTgt_raw, fs_t] = audioread(path_Tgt); % read Target file
    [wavMsk_raw, fs_m] = audioread(path_Msk); % read Masker file
    
    dur_Tgt = timerange_Tgt(2) - timerange_Tgt(1);
    dur_Msk = timerange_Msk(2) - timerange_Msk(1);

    if fs_t~=fs_m
        resample(dur_Msk, fs_t,fs_m);
    end
    
    if dur_Tgt < dur_Msk
        wavTgt = wavTgt_raw(timerange_Tgt(1)*fs_t+1 : timerange_Tgt(2)*fs_t);
        wavMsk = wavMsk_raw(timerange_Msk(1)*fs_t+1 : (timerange_Msk(1)+dur_Tgt)*fs_t);
    else
        wavTgt = wavTgt_raw(timerange_Tgt(1)*fs_t+1 : timerange_Tgt(2)*fs_t);
        wavMsk = [wavMsk_raw(timerange_Msk(1)*fs_t+1 : timerange_Msk(2)*fs_t), ...
                    wavMsk_raw(timerange_Msk(1)+1 : (timerange_Msk(1)+dur_Tgt - dur_Msk)*fs_t)];
    end
    
    %%% S/N configuration %%%
    SNR = 0;
    L = rms(wavMsk) * 10^(SNR/20)/ rms(wavTgt); %SNR = 20*log10(rms(target)/rms(masker))
    wavTgt = L* wavTgt; 
    
    wavTgt = fadein(d_fifo, wavTgt, fs_t);
    stimulus(:,1) = fadeout(d_fifo, wavTgt, fs_t);
    wavMsk = fadein(d_fifo, wavMsk, fs_t);
    stimulus(:,2) = fadeout(d_fifo, wavMsk, fs_t);
    
    disp('finish stimuli preparation')
end
