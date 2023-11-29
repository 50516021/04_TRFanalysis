%%% mTRF analysis - step 6 mTRF SNR comparison %%% 
%%% - SNRs
%%%
%%% required Add-ons
%%% - 
%%% - 
%%% required functions
%%% - 
%%%
%%% required setting files
%%% - 

%%% v1  
%%% 20230928 TRF signal comparison
%%% 20231030 SNR data save section

clearvars; 
close all;

addpath('../'); %add path above
addpath('../../02_EEGanalysis'); %add path of EEGanalysis

OSflag = OSdetection_v1;

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 
outfolder_mTRFfig_step6 = strcat(outfolder, 'mTRF_fig/step6/');
mkdir(outfolder_mTRFfig_step6)
outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl/');
mkdir(outfolder_mTRFmdl)

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


%% data preparation

stim_tag = ["Target", "Masker", "Mixed"]; 
stim_dur = [1 2 3 4 5 6 7 8 9 10]*60; %duration for analysis
stim_dur_load = 10*60; %loading stimuli duration

%%% chanel info
chs = ["Fz", "Cz"];
% chs = ["Fz"];

%%% make legend array for stimuli
for i = 1:length(stim_dur)
    stim_dur_leg(i) = strcat(string(num2str(stim_dur(i)/60)), "min");
    stim_dur_lbl(i) = string(num2str(stim_dur(i)/60));
end

%% load models

numStimtag = size(stim_tag,2);
numStimdur = length(stim_dur);


for j =1: numStimtag
    for k = 1:length(chs)
        figure;
        sgtitle(sprintf('Speech TRF stimulus: %s,  %s ', stim_tag(j), chs(k)))
        for i = 1:numStimdur

            filename = sprintf('mTRF_%s_%0.0fs_%s', stim_tag(j), stim_dur(i), chs(k));
            filename_mdl = strcat(outfolder_mTRFmdl, 'model_', filename, '.mat');
            load(filename_mdl);
            disp(strcat(filename, ' has been loaded'))
            [x, y] = mTRFplot_pros(model,'trf','all',k,[-50,350]);
            if i < numStimdur/4;       plot(x, y, '-'); hold on;
            elseif i < numStimdur*2/4; plot(x, y, '--'); hold on;
            elseif i < numStimdur*3/4; plot(x, y, ':'); hold on;
            else;                      plot(x, y, '.-'); hold on;
            end
            SNR_TRF(i,k,j) = max(y)/rms(y);
            %index: SNR_TRF(stimuli duration, channel, stimuli tag)
        end    
        colororder(parula(numStimdur)); %color gradation
        legend(stim_dur_leg)
        xlim([-50,350])
        ylim([(-1.5)*10^(-3), (1.5)*10^(-3)])
        xlabel('Time lag (ms)')
        ylabel('Amplitude (a.u.)')
        grid on;
        filename_pdf = strcat(outfolder_mTRFfig_step6, sprintf('TRF_%s_%s', stim_tag(j), chs(k)), '.pdf');
        saveas(gcf, filename_pdf)

    end          
end

%% SNR figures

co  = 0;
figure;

for j =1: numStimtag
    for k = 1:length(chs)
        co = co + 1;
        plot(squeeze(SNR_TRF(:,k,j))); hold on;
        leg_SNR(co) = string(sprintf('%s, %s', stim_tag(j), chs(k)));        
    end 
end

legend(leg_SNR)
xticks(1:numStimdur)
xticklabels(stim_dur_lbl)
xlabel('duration [min]');
ylabel('SNR between RMS and Maximum value of TRF');
grid on;
sgtitle(sprintf('SNR of Speech TRF stimulus %s', experiment_name))
filename_pdf = strcat(outfolder_mTRFfig_step6, sprintf('SNR_%s', experiment_name), '.pdf');
saveas(gcf, filename_pdf)

%% save SNR data 

filename_data = strcat(outfolder, sprintf('step6_SNR_%s', experiment_name), '.mat');
save(filename_data, 'SNR_TRF')

%% mTRF processing for plot function
function [x, y] = mTRFplot_pros(model,type,feat,chan,xlims,avgfeat,avgchan)

% Set default values
if nargin < 2 || isempty(type)
    type = 'trf';
end
if model.Dir == -1
    model.w = permute(model.w,[3,2,1]);
end
if nargin < 3 || isempty(feat) || strcmpi(feat,'all')
    feat = 1:size(model.w,1);
end
if nargin < 4 || isempty(chan) || strcmpi(chan,'all')
    chan = 1:size(model.w,3);
end
if nargin < 5 || isempty(xlims)
    xlims = [model.t(1),model.t(end)];
end
if nargin < 6 || isempty(avgfeat)
    avgfeat = 1;
end
if nargin < 7 || isempty(avgchan)
    avgchan = 0;
end

% Define lags
switch type
    case {'mtrf','mgfp'}
        [~,idx1] = min(abs(model.t-xlims(1)));
        [~,idx2] = min(abs(model.t-xlims(2)));
        lags = idx1:idx2;
end

% Define features and channels
model.w = model.w(feat,:,chan);

% Average features
switch type
    case {'trf','gfp'}
        if avgfeat
            model.w = mean(model.w,1);
        end
end

% Average channels
switch type
    case {'trf','mtrf'}
        if avgchan
            model.w = mean(model.w,3);
        end
    case {'gfp','mgfp'}
        model.w = std(model.w,1,3);
end

x = model.t; 
y = model.w;

end