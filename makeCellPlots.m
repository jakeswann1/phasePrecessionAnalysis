function makeCellPlots(Res, prms)
% Jake's function to plot individual cells with multiple potential place
% fields based on lines 164- and onwards of makePlotsPhaseAnalysis.m (plot
% all cells)

%create axes
offSet   = [20 20];
plotSep  = [20 20];
plotSize = [200 100];
screenSz = get(0,'ScreenSize');
       
hScroll = scanpix.plot.createScrollPlot( [0.1*screenSz(3) 0.1*screenSz(4) 3*plotSize(1)+8*plotSep(1) 0.8*screenSz(4) ] );

ageBins = [Res.age Res.age];
animal  = char(Res.uniqueID);

title(char(Res.uniqueID));

% grab data
mapsAll      = { vertcat(Res.mapsShift(ageInd,1)), vertcat(Res.mapsShift(ageInd,2)) };
allPhaseData = { vertcat(Res.allPhaseData(ageInd,1)), vertcat(Res.allPhaseData(ageInd,2)) };
allCircReg   = vertcat(Res.circReg(ageInd));


