%%% mTRF analysis - step 7 mTRF SNR comparison between different participants v3 %%% 
%%% - evoked responce
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

clearvars; 
close all;

%% parameters

addpath('../'); %add path above

OSflag = OSdetection_v1;
figurepath = "figure/";

stim_tag    = ["Target", "Masker", "Mixed"]; 
numTag      = length(stim_tag);
chs         = ["Fz", "Cz"];
numCh       = length(chs);
instruction = ["Left", "Right"];
numInst     = length(instruction);

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
    SNRfile = ls(strcat(foldTemp, 'step6_SNR_s*')); %find responce file
    if OSflag(1) == "1"   
        SNRfile = SNRfile(1:end-1); %extract unnecessary charactar
        SNRdata(i,:,:,:) = struct2array(load(SNRfile)); %SNR data   
    elseif OSflag(1) == "2"
        SNRdata(i,:,:,:) = struct2array(load([foldTemp SNRfile])); %SNR data
    end
    %index:SNRdata(subject, stimuli duration, channel, stimuli tag)

end

%index: SNRdata(participants, stimuli duration, channel, stimuli tag)

%% data preparation 

SNRdata_mean  =  squeeze(mean(SNRdata, 1)); %mean for all participants 

%%% make the table for box plot
co = 0;
for i = 1:numCh
    for j = 1:numTag
        for l = 1:Snum
            co = co + 1;
            
            Channel(co) = chs(i);
            StimTag(co) = stim_tag(j);
            Vari(co)    = string(sprintf('%s, %s', stim_tag(j), chs(i))); 
            subj(co)    = subList.ID(l);
            SNR(co)     = SNRdata(l,k,i,j) ;
           
        end
    end
end


SNR_tbl = table(StimDur', Channel', StimTag', Vari', subj', SNR');
SNR_tbl.Properties.VariableNames = ["StimDur", "Channel", "StimTag", "Variation", "subject", "SNR"];

%% draw figures 

yrange = [1 5];
cols =  ["#0072BD", "#D95319", "#EDB120"];

for i = 1:numCh
    figure;
        
    SNR_tbl_temp = SNR_tbl(SNR_tbl.Channel == chs(i),:);
    SNR_tbl_temp.StimDur = categorical(SNR_tbl_temp.StimDur,stim_dur_lbl);
    groupdata = reordercats(categorical(SNR_tbl_temp.StimTag), stim_tag);
    boxchart(SNR_tbl_temp.StimDur,SNR_tbl_temp.SNR,'GroupByColor',groupdata); hold on;

    for j = 1:numTag  
        plot(squeeze(SNRdata_mean(:,i,j)),'-o', 'Color',cols(j));    hold on; %average lines
        leg_SNR(j) = string(sprintf('%s, %s', stim_tag(j), chs(i)));  
    end

    sgtitle(strcat('Average of TRF SNR:', chs(i), ' (N=', string(Snum), ')'))
    legend([leg_SNR leg_SNR],'Location','northwest')
    ylim(yrange)
    xlabel('duration [min]');
    ylabel('SNR between RMS and Maximum value of TRF');
    grid on;
    hold off;
    
    filename_pdf = strcat(figurepath, sprintf('SNRave_%s_%s', chs(i), subjectlist), '.pdf');
    saveas(gcf, filename_pdf)   
end


