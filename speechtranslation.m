%%% (function) speechtranslation %%% 
%%% - make CRM task stimulus
%%%
%%% required Add-ons
%%% - makestimulus.m
%%% - 
%%% required functions
%%% - 
%%% required setting files
%%% - 

%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.

function  spSpectr = speechtranslation(target, fs, Spat, starttime, SNR, numSpk, Mskflag)

% Mskflag = 0; %whether include masker or not
% target = "0000000";
% fs = 48000;
% Spat = 2;
% starttime = 0;
% SNR = 10;
% numSpk = 3;

%% make speech sound
[stimulus, duration] = makestimulus_TRF(target, fs, Spat, starttime, SNR, numSpk);
speech=stimulus(:,1)+Mskflag*(stimulus(:,Spat+1));

%% convert audio
freqRange = [100 4000];

numFilters = 128;
filterBank = gammatoneFilterBank('SampleRate',fs,'NumFilters',numFilters,'FrequencyRange',freqRange);

speechFiltered = filterBank(speech);
speechAnalytic = hilbert(speechFiltered);

clear spSpectr;

for i = 1:16
    spSpectr(:,i) = mean(abs(speechAnalytic(:,(i-1)*8+1:i*8)),2);
end

% speechSpectrogram = mean(reshape(abs(speechAnalytic),[],8),2);
