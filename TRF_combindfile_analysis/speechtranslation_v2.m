% (function) speechtranslation v2 %%% 
% - translate any mono sound arrays into a vector containing the speech
% spectrogram for mTRF
%
% required Add-ons
% - 
% - 
% required functions
% - 
% required setting files
% - 
% 
%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.
% Parameter


function  spSpectr = speechtranslation_v2(stimulus, fs)

% Mskflag = 0; %whether include masker or not
% target = "0000000";
% fs = 48000;
% Spat = 2;
% starttime = 0;
% SNR = 10;
% numSpk = 3;

%% extract one chanel from the sound array
speech=stimulus(:,1);

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

% %% speech sound plot
% 
% fs = 48000;
% t = 0: 1/fs : 2.8-1/fs; %sample to second conversion
% plot(t, stimulus(1:2.8*fs,1)'); hold on;
% %     ylim([-yscale yscale]); 
% %     xlim([-baselinedur shorteststs]); 
% 
% l1 = line([0;0], get(gca, 'ylim')); l1.Color = [0 1 0]; hold on;
% l2 = line([1;1], get(gca, 'ylim')); l2.Color = [0 0 1]; hold on;
% l3 = line([2.1;2.1], get(gca, 'ylim')); l3.Color = [1 0 1];
% title(sprintf('Speech Example'));
% legend([l1, l2, l3], {'Ready', 'color', 'number'});
% xlabel('time[s]');  
% ylabel('Magnitude');  


