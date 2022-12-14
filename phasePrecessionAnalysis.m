% function Res = phasePrecessionAnalysis(varargin)


%% Set ephys parameters
prms.pathMasterSheet  = 'C:\Users\Jake\Documents\data\replayMasterTrialList_reduced.xlsx';
prms.spkWidthThrForCS = 0.3;  % 0.3, 25 and 5 are all good values, from a histogram of the whole dataset, 2018-04-10. Best prediction is spkWidth vs ACMom, after this, rate only removes 2-3 extra cells.
prms.ACMomThrForCS    = 25;   % 25
prms.rateThrForCS     = 5;    % 5
prms.minNSpikesRun    = 50;

prms.thetaBand        = [5 11]; % This sets the limits of the band in which we search for peak theta frequency (following the FFT)

% plotting
prms.ageBins                    = [16 18; 19 21; 22 25; 26 33; 40 40];    
prms.makeSummaryPlots           = 1;
prms.makeSummaryPlotsPerDataSet = 1;
prms.useSigOnly                 = 0;
prms.pType                      = 'circ_reg'; % 'circ_reg', 'circ_corrc'
prms.plotAllCells               = 1;
prms.multipleFieldPlots         = 1;
%params for density maps
prms.sigma                      = [3 2];

% params for thresholding peaks
prms.minPromThreshold           = 0.5;
prms.maxPromThreshold           = 1;
prms.maxRateThreshold           = 0.5;


% % ---------------------------------------------------------------------------------- %
% if ~isempty(varargin)                                                                %
%     if ischar(varargin{1})                                                           %
%         for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
%     elseif isstruct(varargin{1})                                                     %
%         s = varargin{1};   f = fieldnames(s);                                        %
%         for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
%     end                                                                              %
% end                                                                                  %
% % ---------------------------------------------------------------------------------- %

% pre-allocate output table - only first 3 fields
ResPhase = table( "dummy",...
                  inf,...
                  "dummy",...
                  'VariableNames',{'animal','age','uniqueID'} );

%% fetch data to load from spread sheet
[ expInfo ] = scanpix.helpers.readExpInfo( prms.pathMasterSheet, 'exp' );
%[ expInfo ] = scanpix.helpers.readExpInfo( 'D:\data\replayMasterTrialList.xlsx', 'singletrial','sheetN',1,'animal','r572','experiment',2,'trialName','160726i_sqTrack.set' );
c = 1;
for i = 1:length(expInfo.animal)
    
    % this deals with the length of the various tracks - hardcoded atm!
    % (assuming position is scaled at 400ppm) 
    trackInfo = struct('trackLength',cell(length(expInfo.fullPath{i}),1),'trialType',cell(length(expInfo.fullPath{i}),1));
    
    %Square Track
    if ~isempty( regexp(lower(expInfo.Tnames{i}),'sqtrack','once') )
        idx = ~cellfun('isempty',regexp(lower(expInfo.Tnames{i}),'sqtrack','once'));
        [trackInfo(idx).trialType] = deal('sqtrack');
        if expInfo.age{i}(1) < 40
            [trackInfo(idx).trackLength] = deal(62.5); % 62.5cm/arm
        else
            [trackInfo(idx).trackLength]  = deal(100); % 1m/arm
        end
    end
    
    %Linear Track
    if ~isempty( cell2mat( regexp(lower(expInfo.Tnames{i}),'lintrack','once') ) )
        idx = ~cellfun('isempty',regexp(lower(expInfo.Tnames{i}),'lintrack','once'));
        [trackInfo(idx).trialType] = deal('lintrack');
        if expInfo.age{i}(1) < 40
            [trackInfo(idx).trackLength] = deal(175); % 1.75m long track
%             prms.remTrackEnds = 7;
        else
            [trackInfo(idx).trackLength] = deal(240); % 2.4m long track
%             prms.remTrackEnds = 7;
        end
    end
    
    if all(isnan([trackInfo(:).trackLength]))
        warning('No linear track data found for current data set!');
        continue
    end
      
    %% load data
    % first add cut file identifiers to params container...
    prmsMap = scanpix.helpers.defaultParamsContainer('scanpix.dacq');
    prmsMap('cutTag1') = expInfo.cutFile_ID1{i}{1};
    prmsMap('cutTag2') = expInfo.cutFile_ID2{i}{1};
    % now load
    currObj = scanpix.objLoader('dacq',expInfo.fullPath{i}, 'lfp', true, 'objParams',prmsMap );
    [currObj.trialMetaData.trackLength] = trackInfo.trackLength;
    [currObj.trialMetaData.trialType]   = trackInfo.trialType; 
    
    %% find cross-tetrode recorded cells and remove from data
    runInd = find( ~cellfun('isempty', regexp(lower(currObj.trialNames),'track') ), 1, 'first' ); % atm we only look at the first track trial in case there are more than 1
    scanpix.maps.addMaps(currObj,'lin',runInd);
    scanpix.dacqUtils.findCrossRecCells(currObj,runInd, num2cell(currObj.maps.lin{runInd}{1}, 2), 'removeCellsFromData', 1);
    
    %% filter out IN

    spikeProps = scanpix.analysis.getWaveFormProps(currObj); %get waveform properties
    meanScores = nanmean( spikeProps(:,3:end,:), 3); 
    spkWidthInd  = meanScores(:,1) >= prms.spkWidthThrForCS;
    meanRateInd  = meanScores(:,3) <= prms.rateThrForCS;
    acMomentInd  = meanScores(:,2) <= prms.ACMomThrForCS;
    minNSpikeInd = cellfun(@length, currObj.spikeData.spk_Times{runInd})    >= prms.minNSpikesRun;
    isCSCellInd  = spkWidthInd & meanRateInd & acMomentInd & minNSpikeInd;
    scanpix.helpers.deleteData(currObj,'cells',~isCSCellInd);
    
    %% get best EEG
    % Note that we have hard coded here to use either the unfiltered (filter==0) or unfiltered+Notch (filter==1) or lowpass+notch (filter==4) signal 
    eegInd                                = ismember(currObj.lfpData.lfpTet{runInd}, currObj.cell_ID(:,2)) & ( currObj.trialMetaData(runInd).lfp_filter == 0 | currObj.trialMetaData(runInd).lfp_filter == 1 | currObj.trialMetaData(runInd).lfp_filter == 4);
    eegData                               = horzcat(currObj.lfpData.lfp{runInd}{eegInd});
    [peakFreq, maxPower, s2n, bestEEGInd] = scanpix.lfp.lfpPowerSpec(eegData, prms.thetaBand, 'maxFreq', 25,'eegFs',currObj.params('lfpFs'));
    bestEEG                               = eegData(:,bestEEGInd)';
    thetaFr                               = peakFreq( bestEEGInd );
    % filter
    [~, thetaPhase, cycleN]               = scanpix.lfp.lfpFilter(bestEEG, thetaFr);
    eegIn.phase                           = thetaPhase;
    eegIn.cycleN                          = cycleN;
    %% assign spike phase
    % get maps into right format
    temp.rMaps                                                     = currObj.maps.lin{runInd};
    temp.pMaps                                                     = currObj.maps.linPos{runInd};
    % get posdata in right format
    posData.linPos                                                 = currObj.posData.linXY{runInd};
    posData.linDir                                                 = zeros(length(currObj.posData.linXY{runInd}), 1);
    posData.linDir( ~isnan( currObj.posData.linXY{runInd}(:,2) ) ) = 1;
    posData.linDir( ~isnan( currObj.posData.linXY{runInd}(:,3) ) ) = 2;
    posData.speed                                                  = currObj.posData.speed{runInd};
    tempRes = getSpikePhase(currObj.spikeData.spk_Times{runInd}, eegIn, temp, posData,'trackType', currObj.trialMetaData(runInd).trialType,'LFPFs',currObj.params('lfpFs') );
    tempRes.animal   = string(expInfo.animal{i});
    tempRes.age      = expInfo.age{i}(1);
    tempRes.uniqueID = string([expInfo.animal{i} '_P' num2str(expInfo.age{i}(1))]);

       
    if prms.makeSummaryPlotsPerDataSet
      tempDir = fileparts(prms.pathMasterSheet);
       makePlotsPhaseAnalysis(tempRes,'dataset',prms, [tempDir filesep]); 
%        makePlotsPhaseAnalysis(tempRes,'dataset',prms); 
    end
    
    % join tables - not the most elegant, but saves pre-allocating whole
    % table (for now that's okay)    
    ResPhase.animal(c)   = string(expInfo.animal{i});
    ResPhase.age(c)      = expInfo.age{i}(1);
    ResPhase.uniqueID(c) = string([expInfo.animal{i} '_P' num2str(expInfo.age{i}(1))]);
    if c == 1
        ResPhase         = innerjoin(ResPhase(c,:),tempRes,'Keys',{'animal','age','uniqueID'},'LeftVariables',{'animal','age','uniqueID'});
    else
        ResPhase(c,:)    = innerjoin(ResPhase(c,:),tempRes,'Keys',{'animal','age','uniqueID'},'LeftVariables',{'animal','age','uniqueID'});
    end
    c = c + 1;
end

if prms.makeSummaryPlots
    tempDir = fileparts(prms.pathMasterSheet);
    makePlotsPhaseAnalysis(ResPhase,'population',prms, [tempDir filesep]);
end

% end

