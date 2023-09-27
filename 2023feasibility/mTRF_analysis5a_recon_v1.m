%%% mTRF analysis - step 5 stimuli reconstraction %%% 
%%% - reconstract stimiuli from EEG
%%% - based on https://github.com/mickcrosse/mTRF-Toolbox/blob/master/examples/stimulus_reconstruction.m
%%%
%%% required Add-ons
%%% - 
%%% - 
%%% required functions
%%% - makestimulus_mTRF_v1.m
%%% - StimTrans_mTRF_v1.m
%%% - TRFestimation_v1.m
%%%
%%% required setting files
%%% - 

%%% v1  
%%% 20230926 experiment 'experiment_mTRF_feasibility_v2.m'
%%% non neccesarily comes after step5

clearvars; 
close all;

addpath('../'); %add path above
addpath('../../02_EEGanalysis'); %add path of EEGanalysis

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

%%% get filenames
EEGfile = ls([outfolder, 'step4_*']); %find responce file
% EEGfile = EEGfile(1:end-1); %extract unnecessary charactar
load([outfolder EEGfile]); %participant's responces

%%% get meta data (stimuli info)
metadata_file = ls([outfolder '/metadata*']);
% metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
load([outfolder metadata_file]); %participant's responces


%% stimulus preparation

stimulus = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk);
stimulus(:,3) = stimulus(:,1) + stimulus(:,2); %mixed stimulus
stimulidur = timerange_Tgt(2) - timerange_Tgt(1);

%% data preparation

stim_tag = ["Target", "Masker", "Mixed"]; 
% stim_dur = [5*60 10*60 stimulidur]; %stimuli extraction duration [sec] 
stim_dur = [14*60];

%%% stimuli info

fs_Sound = 48000;
fs_Sound_down = 8000; %consider the frequency range of speech 

stimulus_down = resample(stimulus(1:stim_dur*fs_Sound,:),fs_Sound_down,fs_Sound);

for i =1:size(stimulus,2)
    disp('### begin stimuli conversion ###')
    stim(:,:,i) = StimTrans_mTRF_v1(stimulus_down(:,i), fs_Sound_down);
    disp('### conversion done ###')
end

filename_stim = strcat(outfolder_mTRFfig, 'stim_', experiment_name, '.mat');
save(filename_stim,'stim');

resp = epochs(1:stim_dur*fs_EEG,:);

%% STIMULUS_RECONSTRUCTION  Stimulus reconstruction example.
%   STIMULUS_RECONSTRUCTION loads an example dataset and trains a neural
%   decoder that reconstructs stimulus features (speech envelope) from 2
%   minutes of 128-channel EEG data as per Crosse et al. (2016).
%
%   Example data is loaded from SPEECH_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 15 Hz
%       'fs'        the sample rate of STIM and RESP (128Hz)
%       'factor'    the BioSemi EEG normalization factor for converting the
%                   TRF to microvolts (524.288mV / 2^24bits)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Load data
% load('data/speech_data.mat','stim','resp','fs');

% Normalize data
stim = sum(stim,2);
resp = resp/std(resp(:));

% Downsample data
fsNew = 64;
stim = resample(stim,fsNew,fs_Sound_down);
resp = resample(resp,fsNew,fs_EEG);
fs = fsNew;

% ---Cross-validation---

% Generate training/test sets
nfold = 6;
testfold = 1;
[stimtrain,resptrain,stimtest_arr(:,1),resptest] = mTRFpartition(stim(:,1),resp,nfold,testfold);
for i = 2:size(stim_tag,2)
    stimtest_arr(:,i) = stim(1:size(stimtest_arr,1),i); %different stimulus mix
end

% Model hyperparameters
tmin = 0;
tmax = 250;
lambdas = 10.^(-6:2:6);

% Run fast cross-validation
cv = mTRFcrossval(stimtrain,resptrain,fs,-1,tmin,tmax,lambdas,'zeropad',0,'fast',1);

% ---Model training---

% Get optimal hyperparameters
[rmax,idx] = max(mean(cv.r));
lambda = lambdas(idx);
nlambda = length(lambdas);

% Train model
model = mTRFtrain(stimtrain,resptrain,fs,-1,tmin,tmax,lambda,'zeropad',0);

%% ---Model testing---

for i = 1:size(stim_tag,2)
    
    stimtest = stimtest_arr(:,i);

    % Test model
    [pred,test] = mTRFpredict(stimtest,resptest,model,'zeropad',0);
    
    % ---Evaluation---
    
    % Plot CV correlation
    figure
    subplot(2,2,1)
    errorbar(1:nlambda,mean(cv.r),std(cv.r)/sqrt(nfold-1),'linewidth',2)
    set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1])
    title('CV Accuracy')
    xlabel('Regularization (1\times10^\lambda)')
    ylabel('Correlation')
    axis square, grid on
    
    % Plot CV error
    subplot(2,2,2)
    errorbar(1:nlambda,mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
    set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1])
    title('CV Error')
    xlabel('Regularization (1\times10^\lambda)')
    ylabel('MSE')
    axis square, grid on
    
    % Plot reconstruction
    subplot(2,2,3)
    plot((1:length(stimtest))/fs,stimtest,'linewidth',2), hold on
    plot((1:length(pred))/fs,pred,'linewidth',2), hold off
    xlim([0,10])
    title('Reconstruction')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    axis square, grid on
    legend('Orig','Pred')
    
    % Plot test correlation
    subplot(2,2,4)
    bar(1,rmax), hold on
    bar(2,test.r), hold off
    set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'})
    title('Model Performance')
    xlabel('Dataset')
    ylabel('Correlation')
    axis square, grid on

    sgtitle(sprintf('stimuli reconstraction stimulus: %s, duration:%0.0f s', stim_tag(i), stim_dur))
    filename = sprintf('mTRFrecon_%s_%0.0fs', stim_tag(i), stim_dur);
    filename_pdf = strcat(outfolder_mTRFfig, filename, '.pdf');
    % filename_mdl = strcat(outfolder_mTRFmdl, 'model_', filename, '.mat');
    saveas(gcf, filename_pdf)
    % save(filename_mdl,'model');
    disp(strcat(filename, ' has been saved'))

end

save(strcat(outfolder_mTRFfig, sprintf('mTRFrecon_vars_%0.0fs', stim_dur), '.mat'));

%% stimulus preparation function %%%
function stimulus = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk)

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
    stimulus(:,1) = fadeout(d_fifo, wavTgt, fs_t);
    wavMsk = fadein(d_fifo, wavMsk, fs_t);
    stimulus(:,2) = fadeout(d_fifo, wavMsk, fs_t);
    
    disp('finish stimuli preparation')
end

