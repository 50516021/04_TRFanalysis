%%% TRFanalysis_StimReconst_v1 %%% 
%%% - analyse EEG data by refering wav files for function (for any waves)
%%%
%%% required Add-ons
%%% - mTRF Toolbox
%%% - 
%%% required functions
%%% - speechtranslation_v2.m
%%% required setting files
%%% - 

%%% v1  
%%% 02/20/2023 using proccessed data (after step3 of EEG analysis)
%%% v2 
%%% 03/22/2023 accurate timing syncronization (using behavioral data)


function model = TRFanalysis_StimReconst_v1(stim, resp, fsmTRF)

%%

% Normalize and downsample data
stim = resample(sum(stim,2),64,fsmTRF);
resp = resample(resp/std(resp(:)),64,fsmTRF);
fs = 64;

% Partition data into training/test sets
nfold = 6; testTrial = 1;
[strain,rtrain,stest,rtest] = mTRFpartition(stim,resp,nfold,testTrial);

%% 

% Model hyperparameters
Dir = -1; % direction of causality
tmin = 0; % minimum time lag (ms)
tmax = 250; % maximum time lag (ms)
lambda = 10.^(-6:2:6); % regularization parameters

% Run efficient cross-validation
cv = mTRFcrossval(strain,rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1);

%% 

% Find optimal regularization value
[rmax,idx] = max(mean(cv.r));

% Train model
model = mTRFtrain(strain,rtrain,fs,Dir,tmin,tmax,lambda(idx),'zeropad',0);

% Test model
[pred,test] = mTRFpredict(stest,rtest,model,'zeropad',0);

%% plot figure 

% Plot CV accuracy
figure
subplot(2,2,1), errorbar(1:numel(lambda),mean(cv.r),std(cv.r)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(stest))/fs,stest,'linewidth',2), hold on
plot((1:length(pred))/fs,pred,'linewidth',2), hold off, xlim([0,10]), axis square, grid on
title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,rmax), hold on, bar(2,test.r), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
