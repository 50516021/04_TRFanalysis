%%% (function) StimTrans_mTRF_v1 %%% 
%%% - Sound stimuli translation for mTRF analysis
%%%
%%% required Add-ons
%%% - 
%%% required functions
%%% - 
%%% required setting files
%%% - 

%%% v1
%%% 09/20/2023 sound translation for mTRF analysis
%%% 09/23/2023 devide proccesing for long sound data

%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.

function  stim = StimTrans_mTRF_v1(sound, fs)

% fs = 48000; %sampling rate 
% sound = rand(fs*5, 1); %sound matrix

%% convert audio
freqRange = [100 4000];

numFilters = 128;
filterBank = gammatoneFilterBank('SampleRate',fs,'NumFilters',numFilters,'FrequencyRange',freqRange);

if length(sound) > 3e+06 %long sound data 
    chunk_size = (3e+06)/50;
    chunk_co = fix(length(sound)/chunk_size); %chunk count
    soundAnalytic = [];

    for i = 1:chunk_co
        soundFiltered = filterBank(sound(chunk_size*(i-1)+1:chunk_size*i));
        soundAnalytic = [soundAnalytic; hilbert(soundFiltered)];
        fprintf('filtering %d/%d done \n', i, chunk_co)
    end
    if chunk_size*chunk_co ~= length(sound)
        soundFiltered = filterBank(sound(chunk_size*i+1:end));
        soundAnalytic = [soundAnalytic; hilbert(soundFiltered)];
    end
else
    soundFiltered = filterBank(sound);
    soundAnalytic = hilbert(soundFiltered);
end

clear stim;

for i = 1:16
    stim(:,i) = mean(abs(soundAnalytic(:,(i-1)*8+1:i*8)),2);
end

% soundSpectrogram = mean(reshape(abs(soundAnalytic),[],8),2);

%% sound plot

% t = 0: 1/fs : size(sound,1)*fs-1/fs; %sample to second conversion
% plot(t, stimulus(1:2.8*fs,1)'); hold on;
% %     ylim([-yscale yscale]); 
% %     xlim([-baselinedur shorteststs]); 
% 
% title(sprintf('sound Example'));
% xlabel('time[s]');  
% ylabel('Magnitude');  


