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

%% load list

subjectlist = "subjlist_mTRF_ver20231224";
filesubject = strcat('subjectlist/', subjectlist, '.csv');
listname = string(extractBetween(filesubject, '/', '.csv'));
opts = detectImportOptions(filesubject);
subList = readtable(filesubject, opts);
Snum = size(subList,1);

%% folder information

filenames = {};
for i = 1:Snum
    %%%get folder name
    Snum_str = string(subList.ID(i));
    folder = struct2table(dir(strcat('subject/s000', Snum_str, '*')));
    experiment_name = folder.name; %subject (experiment) name
    foldTemp =  string(sprintf('subject/%s/', experiment_name)); %name of the output folder containing the subject's data 
    filenames{i} = foldTemp;

    % load SNR data
    TRFfile = ls(strcat(foldTemp, ['step6_v5_shortTRF_ICA_s000', char(Snum_str), '*'])); %find responce file
    if OSflag(1) == "1"   %MacOS
        TRFfile = TRFfile(1:end-1); %extract unnecessary charactar
    elseif OSflag(1) == "2" %Windows
        TRFfile = [foldTemp TRFfile];
    end
    load(TRFfile)
    TRFdata(:,:,:,:,:,i)   = TRFs;
    %index: TRFdata([TRF samples], [windows], [instruction], [stimuli], [channel], [subj])

end

num_window = size(TRFs,2);
TRFs_match(:,:,1,:,:)   = TRFs(:,:,1,1,:,:); 
TRFs_match(:,:,2,:,:)   = TRFs(:,:,2,2,:,:); 
TRFs_unmatch(:,:,1,:,:) = TRFs(:,:,2,1,:,:); 
TRFs_unmatch(:,:,2,:,:) = TRFs(:,:,1,2,:,:); 
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

figure('Position', [100 100 800 600]);

co = 0;
for k = 1:Num_Ch_ref
    for j =1: numTag-1 %except mix
        for l = 1:num_window

            co = co+1;
            % subplot(2,2,co)
            nexttile;
            %%% Jackknife calculation and plot
            [tscore(k,j), TRF_JK_mt(:,:,l,k,j), TRF_JK_um(:,:,l,k,j), PeakJK_mt(:,l,k,j), PeakindJK_mt(:,l,k,j), PeakJK_um(:,l,k,j), PeakindJK_um(:,l,k,j)] ...
                = jackknife_comp_v1(TRFs_match(:,l,j,k,:), TRFs_unmatch(:,l,j,k,:), fs, x, peakrangeJK(1), peakrangeJK(2)); 
            %index: TRF_JK_mt(TRF samples, subject, stimuli, channel)
    
            % scatter(x(PeakindJK_mt(:,k,j)), PeakJK_mt(:,k,j), "filled"); hold on;
            % scatter(x(PeakindJK_um(:,k,j)), PeakJK_um(:,k,j), "filled"); hold on;
            % yline(0); hold on;
            % xline(0); hold on;
            % 
            % title(strcat(chs(k), " " ,instruction(j), " stimulus, t = ", string(tscore(k,j))))
            % xlim(xrange) 
            % xlabel('Time lag (ms)');
            % ylim(yrange);
            % ylabel('Amplitude (a.u.)')
            % grid on;

        end
    end
end

% legend("matched", "unmatched", "matched std err", "unmatched std err", "peaks of matched", "peaks of unmatched", 'Location','bestoutside')
% lgd = legend;
% lgd.Layout.Tile = 'east';
% sgtitle(sprintf('TRF Jackknifed average (%d sub) and peaks %s', Snum, subjectlist),'interpreter', "latex")
% filename_pdf = strcat(outfolder_mTRFfig_step7, sprintf('TRFJKave_%s', subjectlist), '.pdf');
% saveas(gcf, filename_pdf)


%% part 2 - attentional modulation
% DataLeft = rand([n_time_of_interest,n_ch,n_subj]);  % put your data [time of interest,channel,subj]
% DataRight = rand([n_time_of_interest,n_ch,n_subj]); % put your data [time of interest,channel,subj]

DataLeft = PeakJK_mt(:,l,k,j);
DataRight = PeakJK_um(:,l,k,j);

topoLeft = squeeze(mean(DataLeft));  % mean in the time of interest 
topoRight = squeeze(mean(DataRight)); % mean in the time of interest

for k=1:n_ch
    [h1(k),p1(k),ci,stat] = ttest(topoLeft(k,:), topoRight(k,:));
    tstat(k) = stat.tstat;
end

figure
set(gcf,'position',[700 605 1000 195])
subplot(1,4,1)
topoplot(nanmean(topoLeft,2),locs,'maplimits',caxis,'whitebk','on')
title(sprintf('Attended Left'))
colorbar
subplot(1,4,2)
topoplot(nanmean(topoRight,2),locs,'maplimits',caxis,'whitebk','on')
title(sprintf('Attended Right'))
colorbar
subplot(1,4,3)
topoplot(tstat,locs,'maplimits',taxis,'whitebk','on')
colorbar
title('T score')
subplot(1,4,4)
topoplot(-log10(p1),locs,'maplimits',paxis,'whitebk','on')
colorbar
title('-log10(p)')
