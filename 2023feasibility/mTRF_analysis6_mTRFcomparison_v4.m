%%% mTRF analysis - step 6 mTRF SNR comparison v4 %%% 
%%% - TRF caricuration for individuals
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
%%% 20231130 nonsub data
%%% v3
%%% 20231207 two instructions 'experiment_mTRF_feasibility_v4.m' (non-duratio/slice)
%%% v4
%%% 20240220 peak comparison
%%% 20240227 ICA option
%%% 20240327 re-referenceing of A1 and A2

clearvars; 
close all;

addpath('../'); %add path above
addpath('../../02_EEGanalysis'); %add path of EEGanalysis

OSflag = OSdetection_v1;

%%% reference %%%
refOpt = ["O1O2" "A1A2"]; %options of onsets
prompt = 'Chose re-referencing cannel:';  % prompt message
[refInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[200 200],'ListString',refOpt); % option selection window
refCh = refOpt(refInd);

%%% option %%%
ICAopt = 1; %use ICA data (1) or not (0)
%%%%%%%%%%%%%%

if ICAopt
    if refInd == 1 %O1 and O2
        namekey = 'step4_plotdatav3_refO1O2*'; %O1 and O2
    elseif refInd == 2
        namekey = 'step4_plotdatav3_refA1A2*'; %A1 and A2
    end  
    nameopt = "_ICA_";
else      
    namekey = 'step4_plotdata_*';   
    nameopt = "_";
end

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 

if ICAopt
    outfolder_mTRFfig_step6 = strcat(outfolder, 'mTRF_fig/step6/');
    outfolder_mTRFmdl       = strcat(outfolder, 'mTRF_mdl/');
else
    outfolder_mTRFfig_step6 = strcat(outfolder, 'mTRF_fig_nonsub/step6/');
    outfolder_mTRFmdl       = strcat(outfolder, 'mTRF_mdl_nonsub/');
end

mkdir(outfolder_mTRFfig_step6)
mkdir(outfolder_mTRFmdl)

if OSflag(1) == "1" %Mac
    %%% get filenames
    EEGfile = ls([outfolder, namekey]); %find responce file
    EEGfile = EEGfile(1:end-1); %extract unnecessary charactar
    load(EEGfile); %participant's responces
    
    %%% get meta data (stimuli info)
    metadata_file = ls([outfolder '/metadata*']);
    metadata_file = metadata_file(1:end-1); %extract unnecessary charactar
    load(metadata_file); %participant's responces

elseif OSflag(1) == "2" %Windows
    %%% get filenames
    EEGfile = ls([outfolder, namekey]); %find responce file
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

for k = 1:numCh
    for j =1: numStimtag-1 %except mix
        for l = 1:length(inst_flg)

            filename = sprintf('mTRF_%s_inst%s', stim_tag(j), instruction(l));
            if refInd == 1 
                filename_mdl = strcat(outfolder_mTRFmdl, 'model', nameopt, filename, '.mat');
            elseif refInd == 2
                filename_mdl = strcat(outfolder_mTRFmdl, 'model', nameopt, refCh, filename, '.mat');
            end
            load(filename_mdl);
            disp(strcat(filename, ' has been loaded'))
            [x, y] = mTRFplot_pros(model,'trf','all',k,[-50,350]);

                if j==l %matched instruction
                    TRF_match(:,k,j) = y;
                else %unmatched
                    TRF_unmatch(:,k,j) = y;
                end

            %index: TRFs(TRF, channel, stimuli tag, instruction)
        end    
    end
end

%% peak ratio figures

xrange = [-50, 350];
yrange = [-7e-4, 10e-4];
peakrange = [0, 350]; 

for k = 1:numCh
    for j =1: numStimtag-1 %except mix
        [max_match(k,j), maxind(k,j)] = max(TRF_match(logical((peakrange(1)<x).*(x<peakrange(2))),k,j)); %only in the case x>0 
        maxind(k,j)                   = maxind(k,j) + sum(x<=peakrange(1)); %x must be adjusted
        max_unmatch(k,j)              = max(TRF_unmatch(maxind(k,j),k,j));  
        MaxRatio_plot(k, j)           = max_match(k,j)-max_unmatch(k,j); %subttraction between matched and unmatched
    end
end

figure;
for k = 1:numCh
    scatter(x(maxind(k,:)), MaxRatio_plot(k,:), "filled"); hold on;
end

xline(0); hold on;
legend(chs)
grid on;
xlim(xrange)
xlabel('Time lag (ms)');
xlabel('stimuli');
ylabel('Subtraction of peak Maximum of TRFs');
grid on;
wholetitle = strcat("Subtraction of TRF peaks: ", experiment_name);
sgtitle(wholetitle,'interpreter', "latex")
filename_pdf = strcat(outfolder_mTRFfig_step6, sprintf('TRFpeaksubq_%s_%s', experiment_name, refCh), '.pdf');
saveas(gcf, filename_pdf)

%% TRF figures

figure;

co = 0;
for k = 1:numCh
    for j =1: numStimtag-1 %except mix
        co = co+1;
        subplot(2,2,co)
        plot(x,TRF_match(:,k,j)); hold on;
        plot(x,TRF_unmatch(:,k,j)); hold on;
        xline(x(maxind(k,j)),'--'); hold on; 
        yline(0); hold on;
        xline(0); hold on;

        title(strcat(chs(k), " " ,instruction(j), " stimulus"))
        xlim(xrange)
        xlabel('Time lag (ms)');
        ylim(yrange);
        ylabel('Amplitude (a.u.)')
        grid on;
    end
end

legend("matched","unmatched", "max for matched", 'Location','best')
sgtitle(sprintf('TRF peaks %s', experiment_name),'interpreter', "latex")
filename_pdf = strcat(outfolder_mTRFfig_step6, sprintf('TRF%s%s_%s',nameopt, experiment_name, refCh), '.pdf');
saveas(gcf, filename_pdf)

%% save SNR data 

filename_data = strcat(outfolder, sprintf('step6_matchTRF%s%s_%s', nameopt, refCh, experiment_name), '.mat');
save(filename_data, 'x', 'TRF_match', 'TRF_unmatch')

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