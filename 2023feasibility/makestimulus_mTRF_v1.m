%%%  makestimulus_mTRF_v1 %%% 
%%% - make stimulus (combination) list, repeat stimulus
%%%     
%%%
%%% #required Add-ons
%%% - 
%%% #required functions
%%% - 
%%% #required setting files
%%%
%%%

%%% v1
%%% 09/14/2023 feasibility test for mTRF analysis (corresponding experiment_mTRF_feasibility_v1.m)


function [stimulus, duration] = makestimulus_mTRF_v1(FstSp, numRep, SpatialPat, SNR, Gap, overwrap, dur)

addpath('../'); %use the same functions as spatial attention test

% FstSp = 1; %["Male", "Female"]: The speaker comes first  
% numRep = 10; %: Repetition number
% SpatialPat = 1; %: Spatial Pattern 1:First-Left/2:First-right
% SNR  = 1; %
% Gap = 1.0; %: Gap between stimuli (sec)
% overwrap = 1.0; % (sec)
% dur = 4.0; % (sec): Duration of one iteration 

%% function variable

d_fifo = 10; %duration of fadein/fadeout (ms)

sound_directory = 'stimulus_mTRF/stimulus20230917'; %'stimlus_mTRF/T-REC-P.501-202005-I!!SOFT-ZST-E/Speech signals/Speech and Noise Signals Clause B/American English_clause_B.3.4/';
if FstSp == 1
    SpeakerFile = ["Male 1.wav", "Female 1.wav"];
else
    SpeakerFile = ["Female 1.wav", "Male 1.wav"];
end

fs = 48000;

%% file reading

    [wavFirst_raw, fs_f] = audioread(sprintf('%s/%s',sound_directory, SpeakerFile(1))); % read Male file
    [wavSecond_raw,fs_s] = audioread(sprintf('%s/%s',sound_directory, SpeakerFile(2))); % read Female file

    %%% resampling
    wavFirst   = resample(wavFirst_raw(1:dur*fs_f),  fs, fs_f);
    wavSecond  = resample(wavSecond_raw(1:dur*fs_s), fs, fs_s);

%% audio stimuli preparation

    %%% S/N configuration %%%
    L = rms(wavSecond) * 10^(SNR/20)/ rms(wavFirst); %SNR = 20*log10(rms(target)/rms(masker))
    wavFirst = L* wavFirst;   

    firstraw  = wavFirst(1:dur*fs);
    secondraw = [zeros(Gap*fs,1)' wavSecond']';
    for i = 1:(numRep-1)
        firstraw  = CrossFader(firstraw,  wavFirst,  overwrap*fs);
        secondraw = CrossFader(secondraw, wavSecond, overwrap*fs);
    end

    %%% S/N configuration %%%
    L = rms(wavSecond) * 10^(SNR/20)/ rms(wavFirst); %SNR = 20*log10(rms(target)/rms(masker))
    wavFirst = L* wavFirst;   
    
%% mix Speach raws
    
    %%% blamk sound matrix
    duration = size(secondraw, 1);
    stimulus = zeros(duration, 2);

    %%% mix
    stimulus(:,SpatialPat)   = fadein(d_fifo, [firstraw' zeros(Gap*fs,1)']', fs);
    stimulus(:,3-SpatialPat) = fadeout(d_fifo, secondraw, fs);
    
end