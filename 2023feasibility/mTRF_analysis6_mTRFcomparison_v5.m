%%% mTRF analysis - step 6 mTRF SNR comparison v5 %%% 
%%% - TRF caricuration for individuals
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
%%% v4
%%% 20240220 peak comparison
%%% 20240227 ICA option
%%% v5
%%% 20240320 sliced short term processing (similer to step5-v6), all cannel
%%% options for topology analysis

clearvars; 
close all;

%%%% options %%%%

Opt_ICA = 1; %use ICA data (1) or not (0)
% Opt_Ch  = 1; %use only Fz and Cz(1) or all (2) -> refers (O1+O2)/2

%%%%%%%%%%%%%%%%%

addpath('../'); %add path above for functions
Locationfile_path = '../../02_EEGanalysis/LocationFiles/';

OSflag = OSdetection_v1;

if Opt_ICA
    namekey = 'step4_plotdatav3_*';
    nameopt = "_ICA_";
else      
    namekey = 'step4_plotdata_*';   
    nameopt = "_";
end

%% load location file
eeglab

locstemp     = readlocs([Locationfile_path, 'DSI-24 Channel Locations w.ced']); %channel configuration file for numCh channels (DSI-24)
ch_list_temp = struct2table(locstemp);
ch_list_orgn = cellstr(ch_list_temp.labels);

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 

if Opt_ICA
    outfolder_mTRFfig = strcat(outfolder, 'mTRF_fig/');
    outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl/');
else
    outfolder_mTRFfig = strcat(outfolder, 'mTRF_fig_nonsub/');
    outfolder_mTRFmdl = strcat(outfolder, 'mTRF_mdl_nonsub/');
end

mkdir(outfolder_mTRFfig)
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

% if Opt_ICA; epochs = saveEp; end %use subtracted (ICAed) epoch


%% parameters

stim_tag = ["Left", "Right"]; 
% stim_dur = [5*60 10*60 stimulidur]; %stimuli extraction duration [sec] 

windowsize = 240;  %sliced window size [sec]
windowgap  = 10; %sliced window gap [sec]

%% stimulus preparation

stimfilename = "stimulus_mTRF/stim_mTRfpilot_v4.mat";
ExpVer = contains(experiment_name,'mTRFpilot_v4'); %experiment version check
ExpStim = exist(stimfilename, 'file');

%%% stimuli info
fs_Sound_down = 8000; %consider the frequency range of speech 
stimulidur = timerange_Tgt(2) - timerange_Tgt(1);

if ExpVer && ExpStim
    disp('### loading stimuli ###')
    stim_str = load(stimfilename);
    stim = stim_str.stim;
    disp('### stimuli have been loaded ###')
else
    [stimulus, fs_Sound] = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk);
    stimulus(:,3) = stimulus(:,1) + stimulus(:,2); %mixed stimulus
    
    %%% resample and extract stimuli
    stimulus_downsmpl = resample(stimulus,fs_Sound_down,fs_Sound);
    
    for i =1:size(stimulus,2)
        disp('### begin stimuli conversion ###')
        stim(:,:,i) = StimTrans_mTRF_v1(stimulus_downsmpl(:,i), fs_Sound_down);
        disp('### conversion done ###')
    end
    
    filename_stim = strcat(outfolder, 'stim_', experiment_name, '.mat');
    save(filename_stim,'stim');
end

%% data preparation

%%% Chanel info
Hotch  = [4 8]; % Fz and Cz (don't change after removing cold channels)
Coldch = [14 15]; % O1 and O2
Num_Ch = size(epochs,2); %number of channels
ch_ind_ref = 1:Num_Ch;
ch_ind_ref(Coldch) = []; % index of channels afetr re-referensing
Num_Ch_ref = length(ch_ind_ref);

ch_list_ref = ch_list_orgn(ch_ind_ref);
epochs_ref = epochs(:,ch_ind_ref,:)-(epochs(:,Coldch(1),:)+epochs(:,Coldch(2),:))/2; %subtraction from center channel

instruction = ["Left", "Right"];
situations = ["matched", "unmatched"];

fs_New = 128;

%% mTRF processing and figure plot

outfolder_mTRFfig_short = strcat(outfolder_mTRFfig, "shortterm/");
mkdir(outfolder_mTRFfig_short)
windowNum = fix((stimulidur-windowsize)/windowgap) + 1; %number of windows
TRFrange = [-50,350];
CFTrange = [100,300]; %timerange of crest factor
yticklabel_values = TRFrange(1):50:TRFrange(2); % y axis parameter

Stim_num = size(stim,3)-1; %except mix

for i = 1:length(inst_flg)
    for j =1:Stim_num
        stim_ext = resample(stim, fs_New, fs_Sound_down);
        resp = resample(epochs(1:stimulidur*fs_EEG,:,i), fs_New, fs_EEG);    

        %%% mTRF estimation (figure)
        figure;
        for k = 1:Num_Ch_ref   
            for l = 1:windowNum
                wndw_strt = ((l-1)*windowgap)*fs_New+1;      %start point of window [sec*fs_New]
                wndw_end  = wndw_strt+windowsize*fs_New-1;    %end   point of window [sec]
                stim_wndw = stim_ext(wndw_strt:wndw_end,:,j); %extracted stimuli 
                resp_wndw = resp(wndw_strt:wndw_end,:);       %extracted EEG responce
                
                %%% model training
                model = TRFestimation_v1(stim_wndw, fs_New, resp_wndw, fs_New, 0, 0);
                if Opt_ICA
                    [x, TRFs(:,l,i,j,k)] = mTRFplot_pros(model,'trf','all',k,TRFrange);
                else
                    [x, TRFs(:,l,i,j,k)] = mTRFplot_pros(model,'trf','all',numch(k),TRFrange);
                end
                %index: TRFs([TRF samples], [windows], [instruction], [stimuli], [channel])
                CFind = CFTrange(1)<=x & CFTrange(2)<=x; %crest factor range 
                CreFac(l,i,j,k) = peak2rms(TRFs(CFind,l,i,j,k));
            end
            
            if  k == 1; clims = [-8 9]*10^(-4);
            else ;      clims = [-7 8]*10^(-4); end %color bar rnage
            % clims = [max(max(max(max(max(TRFs))))) min(min(min(min(min(TRFs)))))];
            
            %%% plot colormap %%%
            subplot(Num_Ch/5,5,k)

            imagesc(squeeze(TRFs(:,:,i,j,k)),clims)
            colorbar

            xlabel("Start Time (s)")
            xlabelsNum = 10; %number of labels
            xvalue = linspace(1,stimulidur-windowsize,xlabelsNum);
            xticks(linspace(1,windowNum,xlabelsNum))
            xticklabels(fix(xvalue))

            ylabel('Time lag (ms)')
            basis = 1; %basis of x memories
            yticks_points = (x(basis+1)-x(basis))^(-1)*(yticklabel_values-x(1));
            yticks(yticks_points)
            yticklabels(yticklabel_values)
            ylim([min(yticks_points) max(yticks_points)])

            title(ch_list_ref{k})

        end   
        % set(gcf,'position',[400 300 700 480])
        sgtitle(sprintf('mTRF stimulus: %s, inst: %s, window size:%2.0fs (gap:%2.1fs) ', stim_tag(j), instruction(i), windowsize, windowgap),'interpreter', "latex")     
        filename = sprintf('mTRF_sliced_wd%dgap%d%s_inst%s',windowsize, windowgap, stim_tag(j), instruction(i));
        filename_pdf = strcat(outfolder_mTRFfig_short, filename, nameopt, '.pdf');
        saveas(gcf, filename_pdf)
        % print(filename_pdf,'-dpdf','-fillpage')
        disp(strcat(filename, ' figure has been saved'))
    end
end

%% save TRF data 

datadir_TRF = strcat(outfolder, "step6_v5_data");
mkdir(datadir_TRF)
filename_data = strcat(datadir_TRF, sprintf('step6_v5_shortTRF%s%s_d%dgap%d', nameopt, experiment_name, windowsize, windowgap), '.mat');
save(filename_data, 'x', 'TRFs', 'CreFac', 'ch_ind_ref', 'Num_Ch_ref', 'ch_list_ref')

%% stimulus preparation function %%%
function [stimulus, fs_t] = stimGen(path_Tgt, timerange_Tgt, path_Msk, timerange_Msk)
    %%% 2023213 resample section

    d_fifo = 10; %duration of fadein/fadeout (ms)
    
    disp('preparing stimuli')
    [wavTgt_raw, fs_t] = audioread(path_Tgt); % read Target file
    [wavMsk_raw, fs_m] = audioread(path_Msk); % read Masker file
    
    dur_Tgt = timerange_Tgt(2) - timerange_Tgt(1);
    dur_Msk = timerange_Msk(2) - timerange_Msk(1);

    if fs_t~=fs_m
        resample(dur_Msk, fs_t,fs_m);
    end
    
    if dur_Tgt < dur_Msk
        wavTgt = wavTgt_raw(timerange_Tgt(1)*fs_t+1 : timerange_Tgt(2)*fs_t);
        wavMsk = wavMsk_raw(timerange_Msk(1)*fs_t+1 : (timerange_Msk(1)+dur_Tgt)*fs_t);
    else
        wavTgt = wavTgt_raw(timerange_Tgt(1)*fs_t+1 : timerange_Tgt(2)*fs_t);
        wavMsk = [wavMsk_raw(timerange_Msk(1)*fs_t+1 : timerange_Msk(2)*fs_t), ...
                    wavMsk_raw(timerange_Msk(1)+1 : (timerange_Msk(1)+dur_Tgt - dur_Msk)*fs_t)];
    end
    
    %%% S/N configuration %%%
    SNR = 0;
    L = rms(wavMsk) * 10^(SNR/20)/ rms(wavTgt); %SNR = 20*log10(rms(target)/rms(masker))
    wavTgt = L* wavTgt; 
    
    wavTgt = fadein(d_fifo, wavTgt, fs_t);
    stimulus(:,1) = fadeout(d_fifo, wavTgt, fs_t);
    wavMsk = fadein(d_fifo, wavMsk, fs_t);
    stimulus(:,2) = fadeout(d_fifo, wavMsk, fs_t);
    
    disp('finish stimuli preparation')
end

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