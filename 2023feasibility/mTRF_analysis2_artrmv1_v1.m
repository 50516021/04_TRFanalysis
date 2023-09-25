%%% mTRF EEG analysis for BDF/EDF/CSV - step 2 %%% 
%%% - artifact remove 1
%%% chose masker or target based mat file
%%%
%%% required Add-ons
%%% - EEGLAB
%%% required functions
%%% - 
%%% required setting files
%%% - BiosemiSettingfiles/BioSemiElecCoor_16.txt
%%% - BiosemiSettingfiles/BioSemiElecCoor_32.txt
%%% - BiosemiSettingfiles/BioSemiElecCoor_64.txt
%%% - BiosemiSettingfiles/BioSemiElecCoor_128.txt
%%% - LocationFiles/DSI-24 Channel Locations w.ced

%%% v1  
%%% 20230920 for mTRF (mTRF_analysis1_epoching_v1.m)

clearvars; close all;

addpath('../../02_EEGanalysis'); %add path of EEGanalysis 
addpath('../../02_EEGanalysis/LocationFiles'); %add path of locationfiles

%% parameters
%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 800],'ListString',folders.name); % option selection window
experiment_name = folders.name(foldInd,:); %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name{1}); %name of the output folder containing the subject's data 

%%% get filenames
data_name = experiment_name{1};
fname = strcat(outfolder, 'step1_', data_name, '.mat'); %epoched EEG data file name with its path
disp(['----- Processing: ' char(data_name), ' -----']) %make sure the processing data

%%%%% parsing excluded channels from the subject database (xls file)
% for k=1:size(xlsStrings,1)
%     if strcmp(xlsStrings{k,1},subjID)
%         ChsExcludedStr = xlsStrings{k,23}
%         break
%     end
% end
% 
% if strcmpi(ChsExcludedStr,'None')
    chs = [];
% else
%     chs = parseChNum(ChsExcludedStr)
% end

%% extract noisy trials and channels

load(fname);
% epochs(:,chs,:) = zeros(size(epochs,1),length(chs),size(epochs,3));

mx = squeeze(max(max(abs(epochs))));
close all
figure
set(gcf,'position',[1   216   560   420])
histogram(mx,100)
title(sprintf('%s, trial selection',data_name));
drawnow();
% thresT = str2double(input('thres? ','s'));
% GoodTrials = find(mx<thresT); %good trial index
title(sprintf('%s trials',data_name));

pdfname = append(outfolder, 'step2_', data_name, '_', 'thTr', '.pdf');
print(pdfname,'-dpdf');

mxC = squeeze(max(mean(abs(epochs),3)));
close all
figure
set(gcf,'position',[1   216   560   420])
histogram(mxC,50)
title(sprintf('%s, channel selection',data_name));
drawnow();
% thresC = str2double(input('thres? ','s'));
% BadChannels = find(mxC > thresC); %bad channel index
title(sprintf('%s channels',data_name));

pdfname = append(outfolder, 'step2_', data_name, '.pdf');
print(pdfname,'-dpdf');

%         save(sprintf('epochs_CCT_OT%04i.mat',sessionIDs(subj)),'GoodTrials','BadChannels','chs_removed','HitOrMiss','-append')
save(strcat(outfolder, "step2_", data_name, "_extraction.mat"),'chs')
% save(strcat(outfolder, "step2_", data_name, "_RejectionThresholds_trials_channels"),  'thresT', 'thresC')

%     fs = 256;
%     norder = 256;
%     cf1 = 1; cf2 = 50;
%     epochs = BPFtd(epochs,fs,norder,cf1,cf2);

%% ICA with EEGLAB
eeglab %load EEGLAB

eeg = permute(epochs,[2,1,3]);

%read location files
numCh   = size(epochs,2); %number of chanels
if numCh == 20
    locstemp    = readlocs('LocationFiles/DSI-24 Channel Locations w.ced'); %channel configuration file for numCh channels (DSI-24)
    locstable = struct2table(locstemp); %swap X and Y
    temp = locstable.X;
    locstable.X = locstable.Y;
    locstable.Y = temp;
    locstable.theta = locstable.theta + 90;
    locs=table2struct(locstable);
else
    locfile = strcat('LocationFiles/BioSemiElecCoor_', num2str(numCh), '.txt'); %channel configuration file for numCh channels (Biosemi)
    locs    = readlocs(locfile,'filetype','xyz'); %load channel configuration file 
end
% pop_chanedit(locs); %if needed

EEG = pop_importdata('dataformat','array', ...
    'nbchan',0,'data','eeg','setname','epochs', ...
    'srate',256,'subject',data_name,'pnts',0, ...
    'xmin',-0.5,'chanlocs','locs');
EEG = pop_select( EEG,'nochannel',chs); %remove specific channels
EEG.chs_removed = chs;

%% save files
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
EEG = pop_saveset(EEG, 'filename', [strcat(outfolder, 'step2_', char(data_name), '_afterRejections_ica.set')]);

disp(['----- Processed: ' char(data_name) ' -----']) %make sure the processing data

close all;
