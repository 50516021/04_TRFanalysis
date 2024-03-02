%%% mTRF analysis - step 6 mTRF SNR comparison v3 %%% 
%%% - SNRs
%%%
%%% required Add-ons
%%% - 
%%% - 
%%% required functions
%%% - 
%%%thi
%%% required setting files
%%% - 

%%% v1  
%%% 20230928 TRF signal comparison
%%% 20231030 SNR data save section
%%% 20231130 nonsub data
%%% v3
%%% 20231207 two instructions 'experiment_mTRF_feasibility_v4.m' (non-duratio/slice)

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
outfolder_mTRFfig_step6 = strcat(outfolder, 'mTRF_fig_nonsub/step6/');
mkdir(outfolder_mTRFfig_step6)
% outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl/');
outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl_nonsub/');

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

stim_tag = ["Left", "Right", "Mixed"]; 
chs = ["Fz", "Cz"];
instruction = ["Left", "Right"];

%% load models

numStimtag = size(stim_tag,2);
numCh = length(chs);

for l = 1:length(inst_flg)
    for j =1: numStimtag
        for k = 1:numCh
            sgtitle(sprintf('Speech TRF stimulus: %s,  %s ', stim_tag(j), chs(k)))

            filename = sprintf('mTRF_%s_inst%s', stim_tag(j), instruction(l));
            filename_mdl = strcat(outfolder_mTRFmdl, 'model_', filename, '.mat');
            load(filename_mdl);
            disp(strcat(filename, ' has been loaded'))
            [x, y] = mTRFplot_pros(model,'trf','all',k,[-50,350]);

            SNR_TRF(k,j,l) = max(y)/rms(y);
            %index: SNR_TRF(channel, stimuli tag, instruction)
        end    
    end
end

%% SNR figures

for l = 1:length(inst_flg)
    co  = 0;
    for k = 1:numCh
        for j =1: numStimtag
            co = co + 1;
            SNR_plot(co,l) = squeeze(SNR_TRF(k,j,l));
            x_label(co,l)  = string(sprintf('%s %s', stim_tag(j), chs(k)));
        end
    end
end

figure;
for l = 1:length(inst_flg)
    plot(SNR_plot(:,l)); hold on;
end

legend(strcat(instruction, " instruction"))
xticks(1:numStimtag*numCh)
xticklabels(x_label(:,1))
xlabel('stimuli/channel');
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