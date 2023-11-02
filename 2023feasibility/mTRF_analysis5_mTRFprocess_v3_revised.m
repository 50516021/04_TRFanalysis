%%% mTRF analysis - step 5 mTRF model generation for various stimuli length %%% 
%%% - 
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
%%% 20231011 changed stim file folder to 'outfolder'
%%% 20231014 load common stim file for mTRFpilot_v2 experiments (fs fixed)

clearvars; 
close all;

addpath('../'); %add path above
addpath('../../02_EEGanalysis'); %add path of EEGanalysis

OSflag = OSdetection_v1;
%%%% Be Careful about sampling rate of sound atimuli

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 
outfolder_mTRFfig = strcat(outfolder, 'mTRF_fig/');
mkdir(outfolder_mTRFfig)
outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl/');
mkdir(outfolder_mTRFmdl)

fs_Sound = 48000;

if OSflag(1) == "1"
    %%% get filenames
    EEGfile = ls([outfolder, 'step4_*']); %find responce file
    EEGfile = EEGfile(1:end-1); %extract unnecessary charactar
    load(EEGfile); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
    load(metadata_file); %participant's responces

elseif OSflag(1) == "2"
    %%% get filenames
    EEGfile = ls([outfolder, 'step4_*']); %find responce file
    load([outfolder EEGfile]); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    load([outfolder metadata_file]); %participant's responces
end

%% stimulus preparation

ExpVer = contains(experiment_name,'mTRFpilot_v2'); %experiment version check

%%% stimuli info
fs_Sound_down = 8000; %consider the frequency range of speech 

if ExpVer
    stim_str = load('stimulus_mTRF/stim_mTRfpilot_v2.mat');
    stim = stim_str.stim;
    disp('### stimuli have been loaded ###')
else
    stimulus = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk, fs_Sound);
    stimulus(:,3) = stimulus(:,1) + stimulus(:,2); %mixed stimulus
    stimulidur = timerange_Tgt(2) - timerange_Tgt(1);
    
    %%% resample and extract stimuli
    stimulus_down = resample(stimulus(1:stim_dur_load*fs_Sound,:),fs_Sound_down,fs_Sound);
    
    for i =1:size(stimulus,2)
        disp('### begin stimuli conversion ###')
        stim(:,:,i) = StimTrans_mTRF_v1(stimulus_down(:,i), fs_Sound_down);
        disp('### conversion done ###')
    end
    
    filename_stim = strcat(outfolder, 'stim_', experiment_name, '.mat');
    save(filename_stim,'stim');
end

%% data preparation

stim_tag = ["Target", "Masker", "Mixed"]; 
% stim_dur = [5*60 10*60 stimulidur]; %stimuli extraction duration [sec] 
stim_dur = [1 2 3 4 5 6 7 8 9 10]*60; %duration for analysis
stim_dur_load = 10*60; %loading stimuli duration

%%% chanel info
chs = ["Fz", "Cz"];
% chs = ["Fz"];

resp = epochs(1:stim_dur*fs_EEG,:);

fs_New = 128;

%% mTRf processing

for i = 1:length(stim_dur)
    stim_ext = resample(stim(1:stim_dur(i)*fs_Sound_down,:,:), fs_New, fs_Sound_down);
    EEG = resample(saveEp(1:stim_dur(i)*fs_EEG,:), fs_New, fs_EEG);
    for j =1:size(stim,3)
        
        %%% mTRF estimation
        for k = 1:length(chs)
            model = TRFestimation_v1(stim_ext(:,:,j), fs_New, EEG, fs_New, k);
            
            sgtitle(sprintf('mTRF stimulus: %s, duration:%0.0f s, %s ', stim_tag(j), stim_dur(i), chs(k)))
            filename = sprintf('mTRF_%s_%0.0fs_%s', stim_tag(j), stim_dur(i), chs(k));
            filename_pdf = strcat(outfolder_mTRFfig, filename, '.pdf');
            filename_mdl = strcat(outfolder_mTRFmdl, 'model_', filename, '.mat');
            saveas(gcf, filename_pdf)
            save(filename_mdl,'model');
            disp(strcat(filename, ' has been saved'))

        end    
    end
end

%% stimulus preparation function %%%
function stimulus = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk, fs_New)

    d_fifo = 10; %duration of fadein/fadeout (ms)
    
    disp('preparing stimuli')
    [wavTgt_raw, fs_t] = audioread(path_Tgt); % read Target file
    [wavMsk_raw, fs_m] = audioread(path_Msk); % read Masker file
    
    dur_Tgt = timerange_Tgt(2) - timerange_Tgt(1);
    dur_Msk = timerange_Msk(2) - timerange_Msk(1);
    
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
    stimulus(:,1) = resample(fadeout(d_fifo, wavTgt, fs_t), fs_New, fs_t);
    wavMsk = fadein(d_fifo, wavMsk, fs_t);
    stimulus(:,2) = resample(fadeout(d_fifo, wavMsk, fs_m), fs_New, fs_m);

    disp('finish stimuli preparation')
end
