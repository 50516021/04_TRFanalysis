%%  CrossFader
% - cross-fade for two files with fs sampling rate and L overlapping length
% (reffer https://stackoverflow.com/questions/23012679/how-to-smoothly-connect-two-signals-in-matlab)
% 
% required Add-ons
% - 
% - 
% required functions
% - 
% required setting files
% -
% v1  
% 06/13/2023 for EEG and sound files
% 
% wave1: a waveform forwarding wave 2
% wave2: a waveform following wave1
% L: length of the number of overlapped samples 

function CombWave = CrossFader(wave1, wave2, L)

W = linspace(1,0,L)';                                    

wave1(end-L+1:end) = wave1(end-L+1:end).*W;
wave2(1:L) = wave2(1:L).*(1-W);

CombWave = zeros(size(wave1,1) + size(wave2,1) - L, 1);
CombWave(1:size(wave1,1)) = wave1;

CombWave(end-size(wave2)+1:end) = CombWave(end-size(wave2,1)+1:end) + wave2;

