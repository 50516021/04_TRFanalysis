%%% mTRF analysis - step 7 mTRF SNR comparison between different participants v4 %%% 
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
xrange = [-50, 350];
yrange = [-6e-4, 6e-4];
peakrange = [100, 300]; 

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
    folder = struct2table(dir(strcat('subject/s000', string(subList.ID(i)), '*')));
    experiment_name = folder.name; %subject (experiment) name
    foldTemp =  string(sprintf('subject/%s/', experiment_name)); %name of the output folder containing the subject's data 
    filenames{i} = foldTemp;

    % load SNR data
    TRFfile = ls(strcat(foldTemp, 'step6_matchTRF_s*')); %find responce file
    if OSflag(1) == "1"   %MacOS
        TRFfile = TRFfile(1:end-1); %extract unnecessary charactar
    elseif OSflag(1) == "2" %Windows
        TRFfile = [foldTemp TRFfile];
    end
    load(TRFfile)
    TRFs_match(i,:,:,:)   = TRF_match; 
    TRFs_unmatch(i,:,:,:) = TRF_unmatch; 
    %index:TRFs_match(subject, TRF, channel, stimuli tag)

end

%index: SNRdata(participants, stimuli duration, channel, stimuli tag)

%% peak ratio figures

%%% peak picking
for i = 1:Snum
    for k = 1:numCh
        for j =1: numTag-1 %except mix
            [max_match(i,k,j), maxind_match(i,k,j)]       = max(TRFs_match(i, logical((peakrange(1)<x).*(x<peakrange(2))),k,j)); %only in the case x>0 
            maxind_match(i,k,j)                           = maxind_match(i,k,j) + sum(x<=peakrange(1)); %x must be adjusted
            value_unmatch(i,k,j)                          = max(TRFs_unmatch(i, maxind_match(i,k,j),k,j));  
            [max_unmatch(i,k,j),  maxind_unmatch(i,k,j)]  = max(TRFs_unmatch(i, logical((peakrange(1)<x).*(x<peakrange(2))),k,j)); %only in the case x>0 
            maxind_unmatch(i,k,j)                         = maxind_unmatch(i,k,j) + sum(x<=peakrange(1)); %x must be adjusted
            MaxRatio_plot(i,k, j)                         = max_match(i,k,j)-value_unmatch(i,k,j); %subttraction between matched and unmatched
        end
    end
end

figure;
for k = 1:numCh
    for j =1:numTag-1 %except mix
        scatter(x(maxind_match(:,k,j)), MaxRatio_plot(:,k,j), "filled"); hold on;
        legends(k,j) = strcat(chs(k), " ", instruction(j), "stim");
    end
end

xline(0); hold on;
legend(legends)
grid on;
xlim(xrange)
xlabel('Time lag (ms)');
xlabel('stimuli');
ylabel('Subtraction between peak Maximum of TRF');
grid on;
wholetitle = strcat("Subtraction of TRF peaks: ", subjectlist);
sgtitle(wholetitle,'interpreter', "latex")
filename_pdf = strcat(outfolder_mTRFfig_step7, sprintf('TRFpeaksub_%s', subjectlist), '.pdf');
saveas(gcf, filename_pdf)

%% TRF figures

%%% average TRF
for k = 1:numCh
    for j =1: numTag-1 %except mix
        TRFave_match(:,k,j)   = mean(TRFs_match(:,:,k,j));
        TRFave_unmatch(:,k,j) = mean(TRFs_unmatch(:,:,k,j));
    end
end

figure;

co = 0;
for k = 1:numCh
    for j =1: numTag-1 %except mix
        co = co+1;
        subplot(2,2,co)
        plot(x,TRFave_match(:,k,j));   hold on;
        plot(x,TRFave_unmatch(:,k,j)); hold on;
        scatter(x(maxind_match(:,k,j)),   max_match(:,k,j),   "filled"); hold on;
        scatter(x(maxind_unmatch(:,k,j)), max_unmatch(:,k,j), "filled"); hold on;
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

legend("matched ave","unmatched ave", "max for matched", "max for unmatched", 'Location','best')
sgtitle(sprintf('TRF average (%d sub) and peaks %s', Snum, subjectlist),'interpreter', "latex")
filename_pdf = strcat(outfolder_mTRFfig_step7, sprintf('TRFave_%s', subjectlist), '.pdf');
saveas(gcf, filename_pdf)
