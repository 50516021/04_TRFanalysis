%%% TRFanalysis(temporal title) %%% 
%%% - analyse EEG data by refering wav files 
%%%
%%% required Add-ons
%%% - mTRF Toolbox
%%% - 
%%% required functions
%%% - 
%%% required setting files
%%% - 

%%% v1  
%%% 02/20/2023

%%%%%%%%%!!check!!%%%%%%%%%%%
%%% name, tgArray, fpassoutfolder

% clearvars; 
close all;
fsEEG=256;
fssound=48000; %sampling rate for the sound
%% parameters
% file names
name = 'test20230118_Hwan'; %subject (experiment) name
datafolder =  ['../02_EEGanalysis/subject/Subj_' name '/']; %name of the folder containing the subject's data 
outfolder =  ['subject/Subj_' name '/']; %name of the folder containing the subject's data 
fnameMsk = strcat(datafolder, 'step3_epochs_Msk_', name, '_ICAprocessedAfterRejections.mat'); %Masker based epoched EEG data file name with its path
fnameTgt = strcat(datafolder, 'step3_epochs_Tgt_', name, '_ICAprocessedAfterRejections.mat'); %Target based epoched EEG data file name with its path
titlename = strrep(name, '_', ' '); %name for figure title (replased '_' to blank)

if ~exist(outfolder, 'dir') %check folder existance
    mkdir(outfolder)
end

%%% stimuli datat

Mskflag = 0; %whether include masker or not
target = "0000000";
Spat = 2;
starttime = 0;
SNR = 10;
numSpk = 3;

% filter settings
fpass = [1.5 20]; %frequency of low/hi cut (Hz)
fsFilt = 230; %order of filtering

% trial information
Cols = ["blue", "red", "white", "green"]; %color variations
Nums = 1:8; %number variations
SNRs = [-12 -18]; %signal noise ratio variations

trial = 1:10; 

%% load data

load(sprintf('%sres.mat',datafolder)); %participant's responces
tgArray= table2array(res(:,12)); 

load(fnameTgt) %load Target based EEG
numChTgt = size(epochs_Gd,2); %number of channels on target

speechSpectrogram = speechtranslation(target, fssound, Spat, starttime, SNR, numSpk, Mskflag); 


% [stimulus, duration] = makestimulus(target, fs, Spat, starttime, SNR, numSpk);

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

%% resample data

stim2=resample(speechSpectrogram(1:fssound*2.8,:), 256, fssound);
factor = 524.288/2^24;
%number of channel
chan = 4;

%% mTRF train

% Normalize data
indx = find(res.StimulusCharactor=='0000000'); %index of particular speech pattern 
% resp = resample(double(epochs_Gd(79:end,:,indx(1))),128,256);
resp = double(epochs_Gd(78:end,:,indx(3)));

% Model hyperparameters
tmin = -100;
tmax = 400;
lambda = 0.5;

% Compute model weights
% Note, ridge regression is used instead of Tikhonov regularization to
% avoid cross-channel leakage of the multivariate input features
%%% mTRFtrain(STIM,RESP,FS,DIR,TMIN,TMAX,LAMBDA) %%%
model = mTRFtrain(stim2,resp*factor,256,1,tmin,tmax,lambda,'method','ridge',...
    'split',5,'zeropad',0);

% H = mTRFplot(MODEL,TYPE,FEAT,CHAN,XLIMS)

%% Plot STRF
figure
subplot(2,2,1), mTRFplot(model,'mtrf','all',chan,[-50,350]);
title('Speech STRF (Fz)'), ylabel('Frequency band'), xlabel('')

% Plot GFP
subplot(2,2,2), mTRFplot(model,'mgfp','all','all',[-50,350]);
title('Global Field Power'), xlabel('')

% Plot TRF
subplot(2,2,3), mTRFplot(model,'trf','all',chan,[-50,350]);
title('Speech TRF (Fz)'), ylabel('Amplitude (a.u.)')

% Plot GFP
subplot(2,2,4), mTRFplot(model,'gfp','all','all',[-50,350]);
title('Global Field Power')

sgtitle('Speech: Ready Blue One')
saveas(gcf, strcat(outfolder, 'TRFmodel_', name, '.pdf'))

% 
% %%%%%% stimulus reconstruction %%%%%%%%%%%5
% % Normalize and downsample data
% stim = resample(sum(stim,2),64,fs);
% resp = resample(resp/std(resp(:)),64,fs);
% fs = 64;
% 
% % Partition data into training/test sets
% nfold = 6; testTrial = 1;
% [strain,rtrain,stest,rtest] = mTRFpartition(stim,resp,nfold,testTrial);
% 
% 
% 
% % Model hyperparameters
% Dir = -1; % direction of causality
% tmin = 0; % minimum time lag (ms)
% tmax = 250; % maximum time lag (ms)
% lambda = 10.^(-6:2:6); % regularization parameters
% lambdas = 10.^(-6:2:6);
% nlambda = length(lambdas);
% 
% % Run efficient cross-validation
% cv = mTRFcrossval(strain,rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1);
% 
% 
% 
% % Find optimal regularization value
% [rmax,idx] = max(mean(cv.r));
% 
% % Train model
% model = mTRFtrain(strain,rtrain,fs,Dir,tmin,tmax,lambda(idx),'zeropad',0);
% 
% % Test model
% [pred,test] = mTRFpredict(stest,rtest,model,'zeropad',0);
% 
% 
% % Plot CV accuracy
% figure
% subplot(2,2,1), errorbar(1:numel(lambda),mean(cv.r),std(cv.r)/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')
% 
% % Plot CV error
% subplot(2,2,2), errorbar(1:numel(lambda),mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')
% 
% % Plot reconstruction
% subplot(2,2,3), plot((1:length(stest))/fs,stest,'linewidth',2), hold on
% plot((1:length(pred))/fs,pred,'linewidth',2), hold off, xlim([0,0.5]), axis square, grid on
% title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')
% 
% % Plot test accuracy
% subplot(2,2,4), bar(1,rmax), hold on, bar(2,test.r), hold off
% set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
% title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
% 
% 
% sgtitle('Speech: Ready Blue One')
% saveas(gcf, strcat(outfolder, 'TRF_reconstruction_', name, '.pdf'))
