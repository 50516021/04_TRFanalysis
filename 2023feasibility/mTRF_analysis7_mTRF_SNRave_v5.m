%%% mTRF analysis - step 7 mTRF SNR comparison between different participants v5 %%% 
%%% - TRF calcuration for multiple people
%%%
%%% required Add-ons
%%% - 
%%% - 
%%% required functions
%%% - 
%%%
%%% required setting files
%%% - subjectlist/subjlist_mTRF_*.csv

%%% v1  
%%% 20231030 inbetween subject comparison (for step6-v1)
%%% v3
%%% 20231224 two instructions 'experiment_mTRF_feasibility_v4.m' (non-duratio/slice)
%%% v4 
%%% 20240221 subtraction plot from step6 v4
%%% 20240313 Jacknife
%%% v5
%%% 20240320 topology analysis
%%% 20240417 gif generation

clearvars; 
close all;

%% parameters

addpath('../'); %add path above

OSflag = OSdetection_v1;
figurepath = "figure/";

stim_tag = ["Left", "Right", "Mixed"];
numTag      = length(stim_tag);
chs         = ["Fz", "Cz"];
numCh       = length(chs);
instruction = ["Left", "Right"];
numInst     = length(instruction);

%%% for peak picking
xrange      = [-50, 350];
yrange      = [-6e-4, 9e-4];
peakrange   = [100, 300]; 
peakrangeJK = peakrange; %peak range for Jackknife

outfolder_mTRFfig_step7 = strcat(figurepath, 'step7/');
mkdir(outfolder_mTRFfig_step7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
windowsize = 240;  %sliced window size [sec]
windowgap  = 10; %sliced window gap [sec]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% options

%%% reference %%%
refOpt = ["O1O2" "A1A2"]; %options of onsets
prompt = 'Chose re-referencing cannel:';  % prompt message
[refInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[200 200],'ListString',refOpt); % option selection window
refCh = refOpt(refInd);

proOpt = ["FzCz" "Allch"]; %options of procesing data
prompt = 'Process:';  % prompt message
[proInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[200 200],'ListString',proOpt); % option selection window
proCh = proOpt(proInd);

%%% option %%%
ICAopt = 1; %use ICA data (1) or not (0)
%%%%%%%%%%%%%%

if ICAopt
    if refInd == 1 %O1 and O2
        namekey = 'step6_v5_shortTRF_ICA_refO1O2_'; %O1 and O2
        Coldch = [14 15]; % O1 and O2
    elseif refInd == 2
        namekey = 'step6_v5_shortTRF_ICA_refA1A2_'; %A1 and A2
        Coldch = [9 18]; %A1 and A2
    end  
    nameopt = "_ICA_";
else      
    namekey = 'step6_matchTRF_s*';   
    nameopt = "_";
end

%% chanel info

if proInd == 1 %use subtracted (ICAed) epoch
    chs = ["Fz", "Cz"];
elseif proInd == 2 %use all channelss epoch
    eeglab
    
    locstable  = struct2table(readlocs('../LocationFiles/DSI-24 Channel Locations w.ced')); %channel configuration file for numCh channels (DSI-24)
    chs_temp  = string(locstable.labels);
    Num_Ch_temp = length(chs_temp); %number of channels
    ch_ind_ref = 1:Num_Ch_temp;
    ch_ind_ref(Coldch) = []; % index of channels afetr re-referensing

    chs = chs_temp(ch_ind_ref);

    temp            = locstable.X;
    locstable.X     = locstable.Y;
    locstable.Y     = temp;
    locstable.theta = locstable.theta + 90;

end

numCh = length(chs);

%% load list

subjectlist = "subjlist_mTRF_ver20240327";
filesubject = strcat('subjectlist/', subjectlist, '.csv');
listname = string(extractBetween(filesubject, '/', '.csv'));
opts = detectImportOptions(filesubject);
subList = readtable(filesubject, opts);
Snum = size(subList,1);

%% folder information

filenames = {};
for i = 1:Snum
    %%%get folder name
    if     subList.ID(i) <  100; subdigit = 's000';    % determine subject's digit
    elseif subList.ID(i) >= 100; subdigit = 's00'; end

    folder = struct2table(dir(strcat('subject/', subdigit, string(subList.ID(i)), '*')));
    experiment_name = folder.name; %subject (experiment) name
    foldTemp =  string(sprintf('subject/%s/', experiment_name)); %name of the output folder containing the subject's data 
    filenames{i} = foldTemp;

    % load SNR data
    TRFfile = ls(sprintf('%s%s%s_d%dgap%d.mat',foldTemp, namekey, experiment_name, windowsize, windowgap)); %find responce file
    if OSflag(1) == "1"   %MacOS
        TRFfile = TRFfile(1:end-1); %extract unnecessary charactar
    elseif OSflag(1) == "2" %Windows
        TRFfile = [foldTemp TRFfile];
    end
    load(TRFfile)
    TRFdata(:,:,:,:,:,i)   = TRFs;
    if exist("TRFs_match")
        TRFs_match_all(:,:,:,:,i)   = TRFs_match;
        TRFs_unmatch_all(:,:,:,:,i) = TRFs_unmatch;
    end
    %index: TRFdata([TRF samples], [windows], [instruction], [stimuli], [channel], [subj])

end

num_window = size(TRFs,2);
if ~exist("TRFs_match")
    TRFs_match_all(:,:,1,:,:)   = TRFdata(:,:,1,1,:,:); 
    TRFs_match_all(:,:,2,:,:)   = TRFdata(:,:,2,2,:,:); 
    TRFs_unmatch_all(:,:,1,:,:) = TRFdata(:,:,2,1,:,:); 
    TRFs_unmatch_all(:,:,2,:,:) = TRFdata(:,:,1,2,:,:); 
end
%index: TRFs_match([TRF samples], [windows], [stimuli], [channel], [subj])

%% load sensor location file

locstemp  = readlocs('../LocationFiles/DSI-24 Channel Locations w.ced'); %channel configuration file for numCh channels (DSI-24)
locstable = struct2table(locstemp); %swap X and Y
temp = locstable.X;
locstable.X = locstable.Y;
locstable.Y = temp;
locstable.theta = locstable.theta + 90; %lotate
locs=table2struct(locstable);
locs = locs(ch_ind_ref); %channels afetr re-referensing

%% Jackknife figures

%%% Jackknife
fs = length(x)/(x(end)-x(1))*1000; %sampling frequency (Hz)
% t_stt = peakrangeJK(1)/1000; %peakrange should be sesond, not ms
% t_end = peakrangeJK(2)/1000;

for k = 1:Num_Ch_ref
    for j =1: numTag-1 %except mix
        for l = 1:num_window

            %%% Jackknife calculation and plot
            [tscore(k,j), TRF_JK_mt(:,l,k,j,:), TRF_JK_um(:,l,k,j,:), PeakJK_mt(l,j,k,:), PeakindJK_mt(l,j,k,:), PeakJK_um(l,j,k,:), PeakindJK_um(l,j,k,:)] ...
                = jackknife_comp_v1(squeeze(TRFs_match_all(:,l,j,k,:)), squeeze(TRFs_unmatch_all(:,l,j,k,:)), fs, x, peakrangeJK(1), peakrangeJK(2)); 
            %index: TRF_JK_mt([TRF samples], [windows], [stimuli], [channel], [subj])
            %index: PeakJK_mt([windows], [stimuli], [channel], [subj])
    
        end
    end
end

%% Topology figures and gif

if proInd == 2

outfolder_gif = strcat(outfolder_mTRFfig_step7, "gif/");
mkdir(outfolder_gif);

for j =1: numTag-1 %except mix
    filename_topo = strcat(outfolder_gif, sprintf('TRFTopoJK_%s_%s_%s_wd%dgap%d', refCh, subjectlist, stim_tag(j), windowsize, windowgap)); 
    filename_topogif = strcat(filename_topo, '.gif');
    for l = 1:num_window
        TRF_Topology_v1(squeeze(PeakJK_mt(l,j,:,:)), squeeze(PeakJK_um(l,j,:,:)), locs)

        window_stt = (l-1)*windowgap;
        window_end = (l-1)*windowgap+windowsize;
    
        sgtitle(sprintf('TRF Topology (JN, %d sub) and peaks %s, stim:%s, window %d - %d s\n', Snum, subjectlist, stim_tag(j), window_stt, window_end'),'interpreter', "latex")
        % filename_pdf = strcat(outfolder_mTRFfig_step7, sprintf('TRFTopoJK_%s_%s_%s', refCh, subjectlist, stim_tag(j)), '.pdf');
        % exportgraphics(gcf,filename_pdf','Resolution',300)
        exportgraphics(gcf, filename_topogif, Append=true);
    end
end

end
