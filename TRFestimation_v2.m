%%% TRFestimation_v2 %%% 
%%% - estimate TRF from arbitrary sound and EEG
%%%
%%% required Add-ons
%%% - mTRF Toolbox
%%% - 
%%% required functions
%%% - 
%%% required setting files
%%% - 

%%% v1  
%%% 20230921 simple function
%%% 20231129 function just for training
%%% v2
%%% 20231213 lambda optimizationz


function model = TRFestimation_v2(stim_original, fs_Sound, EEG, fs_EEG, chan, model, lambda)

%% parameters

% fs_EEG = 256; %sampling rate for EEG (already resampled)
% fs_Sound = 48000; %sampling rate for the sound
fs_mTRF = fs_EEG; %sampling rate for mTRF
% chan = 1; %EEG channel, if chan = 0, just do training

%% resample data
stim=resample(stim_original, fs_mTRF, fs_Sound);
factor = 524.288/2^24; %picked from the developper's instruction
resp = EEG*factor;

%% mTRF train
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

% Model hyperparameters (see bottom)
tmin = -100;
tmax = 400;
% lambda = 0.5;
Dir = 1;

% Compute model weights
% Note, ridge regression is used instead of Tikhonov regularization to
% avoid cross-channel leakage of the multivariate input features
%%% mTRFtrain(STIM,RESP,FS,DIR,TMIN,TMAX,LAMBDA) %%%

if ~isstruct(model)
    model = mTRFtrain(stim,resp,fs_mTRF,Dir,tmin,tmax,lambda,'method','ridge',...
    'split',5,'zeropad',0);
end

%% Plot figure 

if chan~=0

    % Plot STRF
    figure
    subplot(2,2,1), mTRFplot(model,'mtrf','all',chan,[-50,350]);
    title(sprintf('Speech STRF')), ylabel('Frequency band'), xlabel('')
    
    % Plot GFP
    subplot(2,2,2), mTRFplot(model,'mgfp','all','all',[-50,350]);
    title('Global Field Power'), xlabel('')
    
    % Plot TRF
    subplot(2,2,3), mTRFplot(model,'trf','all',chan,[-50,350]);
    title(sprintf('Speech TRF')), ylabel('Amplitude (a.u.)')
    
    % Plot GFP
    subplot(2,2,4), mTRFplot(model,'gfp','all','all',[-50,350]);
    title('Global Field Power')

end

% figtitle = sprintf('Speech:%s-%s-%s, SNR:%d, SpPat:%s, Ch:%s, variation:%s', Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt), Chls{chan}, var);
% sgtitle(figtitle)
% pdfname = sprintf('%s/TRF_%s_Spch-%s-%s-%s_SNR%d_Sp%s_var%s_Ch%s.pdf', outfolderTrial, name, Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt), var, Chls{chan});
% saveas(gcf, pdfname)


%   MTRFTRAIN returns the model in a structure with the following fields:
%       'w'         -- normalized model weights (xvar-by-nlag-by-yvar)
%       'b'         -- normalized bias term (1-by-nlag-by-yvar)
%       't'         -- time lags (ms)
%       'fs'        -- sample rate (Hz)
%       'Dir'       -- direction of causality (forward=1, backward=-1)
%       'type'      -- type of model (multi-lag, single-lag)



