%%% mTRF Topology plot v1 %%% 
%%% - Topology anakysis of TRF
%%%
%%% required Add-ons
%%% - EEGLAB
%%% required functions
%%% - 
%%% required setting files
%%% - 

%%% v1  
%%% 20240320 TRF Topology analysis function (only for the fixed time length, desigened for step7-v4)

function [] = TRF_Topology_v1(TRFpk1, TRFpk2, locs)

%index: TRFpk([channel], [subj]) --- TRF peak files

%% parameters 

caxis = [-2e-4 5e-4];
taxis = [-3 3];
paxis = [0 2];

n_ch   = size(TRFpk1,1); %numbers of Channels
n_subj = size(TRFpk1,2); %numbers of Subjects

%% part 2 - attentional modulation

% index: TRFpk1[channel,subj]

if n_subj>1 %multiple subjects

    for k=1:n_ch
        [h1(k),p1(k),ci,stat] = ttest(TRFpk1(k,:), TRFpk2(k,:));
        tstat(k) = stat.tstat;
    end
    
    TRFpk1 = nanmean(TRFpk1,2);
    TRFpk2 = nanmean(TRFpk2,2);
    
    figure;
    
    set(gcf,'position',[700 605 1000 195])
    subplot(1,4,1)
    topoplot(TRFpk1,locs,'maplimits',caxis,'whitebk','on') 
    title(sprintf('Matched'))
    colorbar
    subplot(1,4,2)
    topoplot(TRFpk2,locs,'maplimits',caxis,'whitebk','on')
    title(sprintf('Unmatched'))
    colorbar
    subplot(1,4,3)
    topoplot(tstat,locs,'maplimits',taxis,'whitebk','on')
    colorbar
    title('T score')
    subplot(1,4,4)
    topoplot(-log10(p1),locs,'maplimits',paxis,'whitebk','on')
    colorbar
    title('-log10(p)')

else %single subject
    
    figure;
    
    set(gcf,'position',[700 605 500 195])
    subplot(1,2,1)
    topoplot(TRFpk1,locs,'maplimits',caxis,'whitebk','on') 
    title(sprintf('Matched'))
    colorbar
    subplot(1,2,2)
    topoplot(TRFpk2,locs,'maplimits',caxis,'whitebk','on')
    title(sprintf('Unmatched'))
    colorbar

end

