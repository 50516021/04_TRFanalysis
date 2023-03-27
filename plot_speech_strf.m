function plot_speech_strf
%PLOT_SPEECH_STRF  Plot example speech STRF.
%   PLOT_SPEECH_STRF loads an example dataset, estimates and plots a speech
%   STRF and the global field power (GFP) from 2 minutes of 128-channel EEG
%   data as per Di Liberto et al. (2015).
%
%   Example data is loaded from SPEECH_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 15 Hertz
%       'fs'        the sample rate of STIM and RESP (128Hz)
%       'factor'    the BioSemi EEG normalization factor for computing the
%                   STRF in microvolts (524.288mV / 2^24bits)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Di Liberto GM, O'Sullivan JA, Lalor EC (2015) Low-Frequency
%          Cortical Entrainment to Speech Reflects Phoneme-Level
%          Processing. Curr Biol 25:1-9.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

close all;
fsEEG=256;
fs=48000;
%% parameters
% file names
name = 'test20230118_Hwan'; %subject (experiment) name
datafolder =  ['../EEGanalysis/subject/Subj_' name '/']; %name of the folder containing the subject's data 
outfolder =  ['subject/Subj_' name '/']; %name of the folder containing the subject's data 
fnameMsk = strcat(datafolder, 'step3_epochs_Msk_', name, '_ICAprocessedAfterRejections.mat'); %Masker based epoched EEG data file name with its path
fnameTgt = strcat(datafolder, 'step3_epochs_Tgt_', name, '_ICAprocessedAfterRejections.mat'); %Target based epoched EEG data file name with its path
titlename = strrep(name, '_', ' '); %name for figure title (replased '_' to blank)

if ~exist(outfolder, 'dir') %check folder existance
    mkdir(outfolder)
end

load(sprintf('%sres.mat',datafolder)); %participant's responces
tgArray= table2array(res(:,12)); %Hwan-112:end, SY-19:160, Minoru-1:113

% filter settings
fpass = [1.5 20]; %frequency of low/hi cut (Hz)
fsFilt = 230; %order of filtering

% trial information
Cols = ["blue", "red", "white", "green"]; %color variations
Nums = 1:8; %number variations
SNRs = [-12 -18]; %signal noise ratio variations
Spat = ["front", "back"]; %spatial pattern variations



stim=resample(speechSpectrogram(1:fs*2.8,:), 256,fs);
factor = 524.288/2^24;
%number of channel
chan = 16;

% Normalize data
indx = find(res.StimulusCharactor=='0000000'); %index of particular speech pattern 
resp = double(epochs_Gd(79:end,:,indx(1)));

% Normalize data
resp = factor*resp;

% Model hyperparameters
tmin = -100;
tmax = 400;
lambda = 0.5;

% Compute model weights
% Note, ridge regression is used instead of Tikhonov regularization to
% avoid cross-channel leakage of the multivariate input features
model = mTRFtrain(stim,resp,fs,1,tmin,tmax,lambda,'method','ridge',...
    'split',5,'zeropad',0);

% Plot STRF
figure, subplot(1,2,1)
mTRFplot(model,'mtrf','all',85,[-50,350]);
title('Speech STRF (Fz)')
ylabel('Frequency band')

% Plot GFP
subplot(1,2,2)
mTRFplot(model,'mgfp','all','all',[-50,350]);
title('Global Field Power')