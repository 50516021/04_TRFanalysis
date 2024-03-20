

caxis = [-1.2 1.2];
taxis = [-3 3];
paxis = [0 2];

n_ch = 20; %?
n_subj = 5; 
n_time_of_interest = 10; % peak in TRF
n_time_segment = 5; % segment of sliding time window

%% part 1 - time course
DataTimeCourse =  rand([n_time_segment,n_ch,n_subj]); % put your data [time windows,channel,subj]
locs  = readlocs('LocationFiles/DSI-24_ChannelLocations.ced'); %channel configuration file for numCh channels (DSI-24)

for k = 1:n_time_segment
    topoTime = squeeze(DataTimeCourse(k,:,:));

    figure
    set(gcf,'position',[700 605 1000 195])
    topoplot(nanmean(topoTime,2),locs,'maplimits',caxis,'whitebk','on')
    title(sprintf('Time frame = %d', k))
    colorbar
    keyboard;
end 

%% part 2 - attentional modulation
DataLeft = rand([n_time_of_interest,n_ch,n_subj]);  % put your data [time of interest,channel,subj]
DataRight = rand([n_time_of_interest,n_ch,n_subj]); % put your data [time of interest,channel,subj]

topoLeft = squeeze(mean(DataLeft));  % mean in the time of interest 
topoRight = squeeze(mean(DataRight)); % mean in the time of interest

for k=1:n_ch
    [h1(k),p1(k),ci,stat] = ttest(topoLeft(k,:), topoRight(k,:));
    tstat(k) = stat.tstat;
end

locs  = readlocs('LocationFiles/DSI-24_ChannelLocations.ced'); %channel configuration file for numCh channels (DSI-24)

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