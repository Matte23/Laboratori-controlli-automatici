function generatePlots(H)
    subplot(1,2,1)
    plotoptions = bodeoptions;
    plotoptions.Grid = 'on';
    %plotoptions.FreqScale = 'linear';
    bodeplot(H,plotoptions)
    subplot(1,2,2)
    [num, den] = tfdata(H, "v");
    nyquist(num,den)
    grid on
end