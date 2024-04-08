

clf;


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