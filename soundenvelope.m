


fs = 48000;

Cols = ["blue", "red", "white", "green"]; %color variations


t = 1: 1/fs : 1.5; %sample to second conversion

for i = 1:4
    subplot(2,2,i)
    plot(t, envelope(spall(fs:fs*1.5,i)))
    title(Cols(i))
    xlabel('time[s]');  

end

sgtitle('Envelope for each Color')