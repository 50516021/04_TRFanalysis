%%% gif_generation v1 %%% 
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
%%% 20240410 gif genartion for short term TRF generated by step6-v5

clearvars;
close all;

%% file info
%example: step6_v5_shortTRF_ICA_refA1A2_s00068_20231121_mTRFpilot_v401_d180gap10.mat

%%%get folder name
folders = struct2table(dir('subject/s*'));
prompt = 'Choose folder name:';  % prompt message
[foldInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',folders.name); % option selection window
experiment_name = folders.name{foldInd,:}; %subject (experiment) name
outfolder =  sprintf('subject/%s/', experiment_name); %name of the output folder containing the subject's data 

%%%choose file to process 
files = struct2table(dir([outfolder, 'step6_v5_shortTRF_*']    ));
prompt = 'Choose file name:';  % prompt message
[fileInd,tf] = listdlg('PromptString',prompt,'SelectionMode','single','ListSize',[400 750],'ListString',files.name); % option selection window
file_path = files.name{fileInd,:}; %subject (experiment) name
file_name = extractAfter(file_path, outfolder);

outfolder_gif = strcat(outfolder, 'mTRF_gif/');
mkdir(outfolder_gif); 

load(file_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data available
% 'TRFs', 'CreFac', 'ch_ind_ref', 'Num_Ch_ref', 'ch_list_ref', 
% 'windowsize', 'windowgap','windowNum',
% 'max_allch_ma','max_allch_um','maxind_ma_ave','maxind_ma_ave',
% 'TRFs_match','TRFs_unmatch'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                               

xrange = [-50, 350];
yrange = [-7e-4, 10e-4];

TRFs_match_Fz = squeeze(TRFs(:,:,1,1,4));

p = plot(nan,nan);
p.XData = x;
dim = [.3 .5 .3 .3];
ano = annotation('textbox',dim,'String',nan,'FitBoxToText','on');

yline(0); hold on;
xline(0); hold on;

% title(strcat(chs(k), " " ,instruction(j), " stimulus"))
xlim(xrange)
xlabel('Time lag (ms)');
ylim(yrange);
ylabel('Amplitude (a.u.)')
grid on;

for n = 1:size(TRFs_match_Fz,2) 
    
    p.YData = TRFs_match_Fz(:,n);
    ano.String = sprintf('n=%d',n);
    exportgraphics(gcf,'testAnimated_1min.gif','Append',true);

end