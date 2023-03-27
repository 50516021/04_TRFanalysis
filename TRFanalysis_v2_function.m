%%% TRFanalysis_v2_function %%% 
%%% - analyse EEG data by refering wav files for function
%%%
%%% required Add-ons
%%% - mTRF Toolbox
%%% - 
%%% required functions
%%% - 
%%% required setting files
%%% - 

%%% v1  
%%% 02/20/2023 using proccessed daat (after step3 of EEG analysis)
%%% v2 
%%% 03/22/2023 accurate timing syncronization (using behavioral data)


function model = TRFanalysis_v2_function(indGen, indCol, indNum, indSNR, indSPt)

clearvars -except indGen indCol indNum indSNR indSPt
close all;

% indGen=1; indCol=1; indNum=1; indSNR=1; indSPt=1;
%% parameters

fsEEG=256; %sampling rate for EEG
fsSound=48000; %sampling rate for the sound
fsmTRF=256; %sampling rate for mTRF

targetdur  = 2.8;  %target duration (sec)
baselinedur= 0.3;  %duration of baseline (sec)

%%% file names
name = 's03_20230317_comptest_Katsu'; %subject (experiment) name
datafolder =  ['../02_EEGanalysis/subject/' name '/']; %name of the folder containing the subject's data 
outfolder =  sprintf('subject/%s/', name); %name of the output folder containing the subject's data 
fnameTgt = strcat(datafolder, 'step3_epochs_Tgt_', name, '_ICAprocessedAfterRejections.mat'); %Target based epoched EEG data file name with its path


if ~exist(outfolder, 'dir') %check folder existance
    mkdir(outfolder)
end

% filter settings
fpass = [1.5 20]; %frequency of low/hi cut (Hz)
fsFilt = 230; %order of filtering

%% trial information

Mskflag = 1; %whether include masker or not
% target = "0000000";
% Spat = 2;
% starttime = 0;
% SNR = 10;
numSpk = 3; %number of speakers (fixed as 3 for the experiment setting)

Gens = ["male", "female"]; %speakers gender variations
Cols = ["blue", "red", "white", "green"]; %color variations
Nums = ["1", "2", "3", "4", "5", "6", "7", "8"]; %number variations
SNRs = ["-12", "-18"]; %signal noise ratio variations
SPts = ["front", "back"]; %spatial pattern variations
% %%% >Gender 
% prompt = 'Choose Speaker Gender:';  % prompt message
% [indGen,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',Gens); % option selection window
% %%% >Color 
% prompt = 'Choose Color:';  % prompt message
% [indCol,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',Cols); % option selection window
% %%% >Number 
% prompt = 'Choose Number:';  % prompt message
% [indNum,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',Nums); % option selection window
% %%% >SNR
% prompt = 'Choose SNR:';  % prompt message
% [indSNR,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',SNRs); % option selection window
SNR = str2double(SNRs(indSNR)); %SNR
%%% >Spatial Pattern
% prompt = 'Choose Spatial Pattern:';  % prompt message
% [indSPt,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',SPts); % option selection window
Spat = indSPt;

%generate taget
target = strcat(string((indGen-1)*4), '000', string(indCol-1), '0', string(indNum-1));

%% load data

%%% get behavioral data
resfile = ls([datafolder '/res*']);
load(resfile(1:end-1)); %participant's responces
tgArray= table2array(res(:,12)); %Hwan-112:end, SY-19:160, Minoru-1:113

%%% EEG data
load(fnameTgt, 'GdTr_final', 'epochs_Gd') %load Target based EEG

%% EEG channel configuration (device detection by the number of channel)
NumCh = size(epochs_Gd,2);

%DSI-24
if NumCh == 20 
    locs    = struct2table(readlocs('../02_EEGanalysis/LocationFiles/DSI-24 Channel Locations w.ced')); %channel configuration file for numCh channels (DSI-24)
%Biosemi
else           
    locfile = strcat('../02_EEGanalysis/LocationFiles/BioSemiElecCoor_', num2str(NumCh), '.txt'); %channel configuration file for numCh channels (Biosemi)
    locs    = struct2table(readlocs(locfile,'filetype','xyz')); %load channel configuration file 
end

Chls    = cellstr(locs.labels); %channel name array

%% determine trial index

% Normalize data
indx = find(res.StimulusCharactor==target); %index of particular target

co = 1; %counter
clear indxTemp;
for i = 1: length(indx)
    if res.SNR(indx(i))==SNR && res.SpatialPosition(indx(i))==Spat %index of particular SNR & Spat           
        indxTemp(co) = indx(i);
        co = co +1;
    end
end

var = 'no'; %variation
if ~any(GdTr_final==indxTemp) %check the existance of the GoodTrial
    var = 'no data'; %variation - there are no eligible trials
    close all; 
    figtitle = sprintf('Speech:%s-%s-%s, SNR:%d, SpPat:%s, variation:%s', Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt), var);
    sgtitle(figtitle)
    pdfname = sprintf('%s/TRF_%s_Spch-%s-%s-%s_SNR%d_Sp%s_var%s.pdf', outfolder, name, Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt), var);
    saveas(gcf, pdfname)
    disp('there are no eligible trials');
    return; %exit
else 
%     prompt = 'Choose Trial Index:';  % prompt message
    co = 1; %counter
    for i=1:size(indxTemp,2)
        if any(GdTr_final==indxTemp(i))
            indxFinalMat(co) = find(GdTr_final==indxTemp(i));
            co = co +1;
        end
    end
     %     [indFind,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[150 100],'ListString',string(indxFinal)); % option selection window
%     indxFinal = indxFinal(indFind);
%     var = char(indFind);
end

for i=1:length(indxFinalMat)
    indxFinal = indxFinalMat(i);
    if co > 2
        var = string(i);
    end

    %%% make trial folder
    foldername = sprintf('Spch-%s-%s-%s_SNR%d_Sp%s_var%s', Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt),var);
    outfolderTrial = [outfolder foldername];
    if ~exist(outfolderTrial, 'dir') %check folder existance
        mkdir(outfolderTrial)
    end
    
    resp = resample(double(epochs_Gd(:,:,indxFinal)),fsmTRF,fsEEG); %extract data and resample
    
    starttime = res.StartTime(indxFinal);

    %%% speech data
    speechSpectrogram = speechtranslation(target, fsSound, Spat, starttime, SNR, numSpk, Mskflag); 
    
    %% resample data
    stim=resample(speechSpectrogram(starttime-baselinedur*fsSound+1:starttime+fsSound*targetdur,:), fsmTRF, fsSound);
    factor = 524.288/2^24;
    %number of channel

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
    
    % Model hyperparameters
    tmin = -100;
    tmax = 400;
    lambda = 0.5;
    
    % Compute model weights
    % Note, ridge regression is used instead of Tikhonov regularization to
    % avoid cross-channel leakage of the multivariate input features
    %%% mTRFtrain(STIM,RESP,FS,DIR,TMIN,TMAX,LAMBDA) %%%
    model = mTRFtrain(stim,resp*factor,256,1,tmin,tmax,lambda,'method','ridge',...
        'split',5,'zeropad',0);
    
    % H = mTRFplot(MODEL,TYPE,FEAT,CHAN,XLIMS)
    for j=1:NumCh
    %     chan = 1; %EEG channel
        chan = j;  
        %% Plot figure 
        
        % Plot STRF
        figure
        subplot(2,2,1), mTRFplot(model,'mtrf','all',chan,[-50,350]);
        title(sprintf('Speech STRF (%s)', Chls{chan})), ylabel('Frequency band'), xlabel('')
        
        % Plot GFP
        subplot(2,2,2), mTRFplot(model,'mgfp','all','all',[-50,350]);
        title('Global Field Power'), xlabel('')
        
        % Plot TRF
        subplot(2,2,3), mTRFplot(model,'trf','all',chan,[-50,350]);
        title(sprintf('Speech TRF (%s)', Chls{chan})), ylabel('Amplitude (a.u.)')
        
        % Plot GFP
        subplot(2,2,4), mTRFplot(model,'gfp','all','all',[-50,350]);
        title('Global Field Power')
    
        figtitle = sprintf('Speech:%s-%s-%s, SNR:%d, SpPat:%s, Ch:%s, variation:%s', Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt), Chls{chan}, var);
        sgtitle(figtitle)
        pdfname = sprintf('%s/TRF_%s_Spch-%s-%s-%s_SNR%d_Sp%s_var%s_Ch%s.pdf', outfolderTrial, name, Gens(indGen), Cols(indCol), Nums(indNum), SNR, SPts(indSPt), var, Chls{chan});
        saveas(gcf, pdfname)
        
    end
end

