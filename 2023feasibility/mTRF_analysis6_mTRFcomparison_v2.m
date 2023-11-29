%%% mTRF analysis - step 6 mTRF SNR comparison %%% 
%%% - time slice
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
%%% v2
%%% 20231101 time slice

clearvars; 
close all;

addpath('../'); %add path above
addpath('../../02_EEGanalysis'); %add path of EEGanalysis

OSflag = OSdetection_v1;

%% folder structure
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 
outfolder_mTRFfig_step6 = strcat(outfolder, 'mTRF_fig/step6/timeslice/');
mkdir(outfolder_mTRFfig_step6)
outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl/');
mkdir(outfolder_mTRFmdl)

% if OSflag(1) == "1"
%     %%% get filenames
%     EEGfile = ls([outfolder, 'step4_*']); %find responce file
%     EEGfile = EEGfile(1:end-1); %extract unnecessary charactar
%     load(EEGfile); %participant's responces
% 
%     %%% get meta data (stimuli info)
%     metadata_file = ls([outfolder '/metadata*']);
%     metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
%     load(metadata_file); %participant's responces
% 
% elseif OSflag(1) == "2"
%     %%% get filenames
%     EEGfile = ls([outfolder, 'step4_*']); %find responce file
%     load([outfolder EEGfile]); %participant's responces
% 
%     %%% get meta data (stimuli info)
%     metadata_file = ls([outfolder '/metadata*']);
%     load([outfolder metadata_file]); %participant's responces
% end

%% parameters
stim_tag = ["Target", "Masker", "Mixed"]; 
stim_dur = [1 2 3 4 5 6 7 8 9]*60; %duration for analysis
stim_dur_load = 10*60; %loading stimuli duration
slice_step = 1*60; %step of slice (sec)
MaxStep = 10*60; %maximum duration of slice (sec)

chs = ["Fz", "Cz"];

%% data preparation

%%% make legend array for stimuli
for i = 1:length(stim_dur)
    stim_dur_leg(i) = strcat(string(num2str(stim_dur(i)/60)), "min");
    stim_dur_lbl(i) = string(num2str(stim_dur(i)/60));
end

%% load models and TRF plot

numStimtag = size(stim_tag,2);
numStimdur = length(stim_dur);
plotPat = ["-", "--", ".-", ":"]; %plot pattern


for j =1: numStimtag
    for k = 1:length(chs)
        for i = 1:numStimdur
            numSlice(i) = (MaxStep-stim_dur(i))/slice_step+1; %number of slice

            figure;
            sgtitle(sprintf('Speech TRF stimulus: %s,  %s - %ds ', stim_tag(j), chs(k), stim_dur(i)/60))
            
            for l = 1:numSlice(i)
                %%% load model
                filename = sprintf('mTRF_%s_%0.0fs_%s_slice%02.f', stim_tag(j), stim_dur(i), chs(k), l);
                filename_mdl = strcat(outfolder_mTRFmdl, 'model_', filename, '.mat');
                load(filename_mdl);
                disp(strcat(filename, ' has been loaded'))

                [x, y] = mTRFplot_pros(model,'trf','all',k,[-50,350]);
                if l < numSlice(i)/4;       plot(x, y, plotPat(1)); hold on;
                elseif l < numSlice(i)*2/4; plot(x, y, plotPat(2)); hold on;
                elseif l < numSlice(i)*3/4; plot(x, y, plotPat(3)); hold on;
                else;                       plot(x, y, plotPat(4)); hold on;
                end
                SNR_TRF(l,i,k,j) = max(y)/rms(y);
                %index: SNR_TRF(slice, stimuli duration, channel, stimuli tag)

                stim_slice_leg(i,l) = strcat("start: ", sprintf('%d sec', (l-1)*60)); %legend of slices (starting point)

            end
            
            colororder(parula(numStimdur)); %color gradation
            legend(stim_slice_leg(i,1:numSlice))
            xlim([-50,350])
            ylim([(-1.5)*10^(-3), (1.5)*10^(-3)])
            xlabel('Time lag (ms)')
            ylabel('Amplitude (a.u.)')
            grid on;
            filename_pdf = strcat(outfolder_mTRFfig_step6, sprintf('TRF_%s_%ds_%s', stim_tag(j), stim_dur(i), chs(k)), '.pdf');
            saveas(gcf, filename_pdf)

        end
    end          
end


%% SNR figures

for i = 1:numStimdur
    co  = 0;
    figure;
    for j =1: numStimtag
        for k = 1:length(chs)
            co = co + 1;
            plot(squeeze(SNR_TRF(1:numSlice(i),i,k,j)), plotPat(j)); hold on;
            leg_SNR(co) = string(sprintf('%s, %s', stim_tag(j), chs(k)));   
        end
    end 

    legend(leg_SNR)
    xticks(1:numSlice(i))
    xticklabels(string(((1:numSlice(i))-1)*60))
    xlabel('start point [sec]');
    ylabel('SNR between RMS and Maximum value of TRF');
    grid on;
    sgtitle(sprintf('SNR of Speech TRF stimulus (%ds) %s', stim_dur(i), experiment_name))
    filename_pdf = strcat(outfolder_mTRFfig_step6, sprintf('SNR_%s_%ds', experiment_name, stim_dur(i)), '.pdf');
    saveas(gcf, filename_pdf)

end

%% save SNR data 

filename_data = strcat(outfolder, sprintf('step6_SNR_timeslice_%s', experiment_name), '.mat');
save(filename_data, 'SNR_TRF')
fprintf('%s has been saved \n', filename_data)

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