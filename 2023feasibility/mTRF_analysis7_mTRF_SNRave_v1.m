%%% mTRF analysis - step 7 mTRF SNR comparison between different participants %%% 
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
%%% 20231030 inbetween ubject comparison (for step6-v1)

clearvars; 
close all;

%% parameters

addpath('../'); %add path above

OSflag = OSdetection_v1;
figurepath = "figure/";

stim_tag = ["Target", "Masker", "Mixed"]; 
numTag = length(stim_tag);
stim_dur = [1 2 3 4 5 6 7 8 9 10]*60; %duration for analysis
numDur = length(stim_dur);
chs = ["Fz", "Cz"];
numCh = length(chs);

%% load list

subjectlist = "subjlist_mTRF_ver20231101";
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

end

%index: SNRdata(participants, stimuli duration, channel, stimuli tag)

%% data preparation 

%%% make legend array for stimuli
for i = 1:length(stim_dur)
    stim_dur_leg(i) = strcat(string(num2str(stim_dur(i)/60)), "min");
    stim_dur_lbl(i) = string(num2str(stim_dur(i)/60));
end

SNRdata_mean  =  squeeze(mean(SNRdata, 1)); %mean for all participants 

%%% make the table for box plot
co = 0;
for k = 1:numDur
    for i = 1:numCh
        for j = 1:numTag
            for l = 1:Snum
                co = co + 1;
                
                StimDur(co) = string(stim_dur_lbl(k));
                Channel(co) = chs(i);
                StimTag(co) = stim_tag(j);
                Vari(co)    = string(sprintf('%s, %s', stim_tag(j), chs(i))); 
                subj(co)    = subList.ID(l);
                SNR(co)     = SNRdata(l,k,i,j) ;
               
            end
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
    boxchart(SNR_tbl_temp.StimDur,SNR_tbl_temp.SNR,'GroupByColor',SNR_tbl_temp.Variation); hold on;
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


