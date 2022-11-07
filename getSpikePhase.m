function [Res] = getSpikePhase(spikeTimes, eegData, maps, posData, varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Params
prms.LFPFs                = 250; % sample rate for eeg 
prms.posFs                = 50;  % sample rate for pos tracking
prms.trackType            = 'sqTrack';
prms.binSize_linMaps      = 10; % in cam pix
prms.phaseDef = 'Skaggs'; % 'Jeewajee'; 'Skaggs'; 'none'
% general cell thresholds 
prms.minSIForMapInclusion = 0;
prms.minPeakRate          = 1;
% thresholds for field detection
prms.fieldDetectMode      = 'allFields'; % 'Pk2LocalMin'; 'MaxArea>Thr'; 'allFields'; pkProminence;
prms.rateThrForField      = 0.4;
prms.sizeThrForField      = 5;
prms.localMinThrForField  = 0.5;

prms.minPromThreshold     = 0.1;
prms.maxPromThreshold     = 1;

% behavioural thresh
prms.speedThrForPhasePrec = 2.5;

% circ regression
prms.mNCyclForRegression  = 2;
prms.minMeth              = 'slope'; % 'slope'; 'slope+off'

% params for density maps
prms.binSzPhaseDeg        = 6;
prms.nPosBinsDensMap      = 20;

% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %

binSzPhase = prms.binSzPhaseDeg * pi/180; % radians

%% assign phases
% 1. assign a phase to each spike (by converting the spike times into an index into the phase vector %

% several possible methods
switch lower(prms.phaseDef)
    case 'skaggs'
        % use phase normalisation from Feng et al. (2015) JNeuro (or rather Skaggs et al. 1996)
        % we can do all data in one go:
        
        allSpikeTimes                            = vertcat(spikeTimes{:});
        allPhaseInd                              = ceil( allSpikeTimes .* prms.LFPFs );
        allPhaseInd( allPhaseInd == 0 )          = 1;
   
        % 1. get start time and length of each theta cycle
        cycleLength                              = accumarray(eegData.cycleN( ~isnan(eegData.cycleN) )',1) ./ prms.LFPFs;  % length in sec of all legit cycles
       
        tempCycleID                              = eegData.cycleN( allPhaseInd ); % cycle ID for each spike
        [cycleID,~,cycleInd]                     = unique(tempCycleID(~isnan(tempCycleID)),'stable'); % list of all cycles + index
        cycleIndByNSamples                       = bsxfun(@eq,eegData.cycleN,cycleID'); % cycleInd x nSamples logical array
        % now we need the start time for each cycle
        [~,startInd]                             = max(cycleIndByNSamples,[],2);
        t0                                       = startInd(cycleInd) ./ prms.LFPFs - 1/prms.LFPFs; % expand and convert to sec (left edges of bins)
        % get the index for each spike
        tempL                                    = cycleLength(cycleID); % expand   

        % 2. interpolate phase 
        allPhase                                 = nan(length(allPhaseInd),1); % need to remember nan's
        allPhase(~isnan(tempCycleID))            = 2*pi .* (allSpikeTimes(~isnan(tempCycleID)) - t0) ./ tempL(cycleInd); % interpolate phase 
        allPhase(allPhase > 2*pi)                = 2*pi;   % see line 64: because we used left edge some samples are >2*pi
        
        % 3. get global phase 0 (phase of max pyramidal cell activity)
        bins                                     = linspace(0,2*pi,361);
        rateHist                                 = histcounts(allPhase,bins,'Normalization','probability' );
        rateHist                                 = imfilter(rateHist,fspecial('gaussian',[1 180],45),'circular'); % smooth with Gaussian (sigma=45deg, so we get very smooth curve)
        globalPhase0                             = mod( bins( ceil( mean( find( rateHist == max(rateHist) ) ) ) ) + 0.5*pi/180, 2*pi );
        allPhase                                 = mod(allPhase - globalPhase0, 2*pi);
        idx                                      = 1; % initialise index for loop - could use accumarray to do the whole shebang in one go
                
    case 'jeewajee'
        % use phase normalisation from Jeewajee et al. (2014) RoySoc - see
        % below for the actual norm.
        k         = fspecial('gaussian',[1 180],45);   % filter for smoothing - as we use 1 deg bin we can keep degrees here
        phaseBins = linspace(0,2*pi,361);     % these are 1 deg bins in radians
end

spikePhase = cell(length(spikeTimes), 1);
for i = 1:length(spikeTimes)
    
    spkIndToPhase                     = ceil( spikeTimes{i} .* prms.LFPFs );
    spkIndToPhase( spkIndToPhase==0 ) = 1;
    
    switch lower(prms.phaseDef)
        case 'skaggs'
            % now just grab data - have done the hard bit already
            spikePhase{i} = allPhase(idx:idx+length(spkIndToPhase)-1);
            idx           = idx+length(spkIndToPhase); % increment index
            
        case 'jeewajee'
            tempBinPhase   = histcounts(eegData.phase( spkIndToPhase ), phaseBins);              % bin phase in to 1 deg bins
            filteredSignal = imfilter(tempBinPhase,k,'circular');                                % smooth with Gaussian (sigma=45deg)
            zeroPhase      = phaseBins(ceil(mean(find(filteredSignal == min(filteredSignal))))); % find minimum (bin with lowest n spikes)
            spikePhase{i}  = mod(eegData.phase( spkIndToPhase ) - zeroPhase,2*pi)';                  % normalise phase, first shift 0...
            % need to normalise mean phase, i.e. set field center phase = 180
            
        case 'none'
            spikePhase{i}  = eegData.phase( spkIndToPhase )'; % do nothing, i.e. raw phase 
     
    end
end

% % 1c. for CurrBiol R1, a reviewer asked about theta phase variability, so here we get the
% %       mean and variance of cycle length, but only for periods when the rat is on a coherent run.
% cycInd          = abs( diff(thetaPhase) ) >= pi;     % This is a logical index for where theta phase wraps.
% cycNumInds      = find( cycInd );
% cycTimes        = [cycNumInds(1:end-1) cycNumInds(2:end)]   ./ LFPSampRate;
% runEpochMask    = Res.speed{RUNTrInd}>=prms.speedThrForPhasePrec  &  Res.linDir{RUNTrInd}~=0;
% runEpochStarts  = (find( diff([0; runEpochMask; 0])==1 ) - 1) ./ posSampRate;
% runEpochEnds    = (find( diff([0; runEpochMask; 0])==-1 ) )   ./ posSampRate;
% winStartInState = bsxfun( @gt, cycTimes(:,1), runEpochStarts' );
% winEndInState   = bsxfun( @lt, cycTimes(:,2), runEpochEnds' );
% winInStateEpoch = winStartInState & winEndInState;          % This is a (nWin, nStateEpoch) logical array, true where a win lies within a state epoch.
% validWinInd     = any(  winInStateEpoch ,   2   );          % However, we are only interested in whether win is valid (not which epoch), so just do ANY along rows.
% cycTimes        = cycTimes( validWinInd, : );
% cycDurs         = diff( cycTimes, 1, 2 ) .* 1000;   % Convert to ms.
% Res.thetaDurMean(RUNTrInd) = mean( cycDurs );
% Res.thetaDurSD(RUNTrInd)   = std( cycDurs );
    

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) **Phase Shift ** For each cell (and treating CW and CCW separately) find the position of the largest field and plot
%     phase versus x for the spikes, while the rat is in the field and running in the appropriate direction.

% grap spatial info for all cells
bits_per_spike = nan(length( spikeTimes ),2); 
for itDir=1:2 
    bits_per_spike(:,itDir) = scanpix.analysis.spatial_info(num2cell(maps.rMaps{itDir+1},2), maps.pMaps(itDir+1,:));
end

% pre-allocate some stuff
maxPixRMap                                                       = size(maps.rMaps{1},2) * prms.binSize_linMaps; % rate map size in cam pix
allPhasePos                                                      = {zeros(0,3), zeros(0,3)};
[m, phi0, rho_CircReg, rho, nSpikesDir, p_CircReg, pVal, pTheta] = deal( nan( length( spikeTimes ), 2 ) ); 
densMapOut                                                       = nan(floor(360/prms.binSzPhaseDeg),prms.nPosBinsDensMap,length(spikeTimes),2);
mapsOut                                                          = {nan(length(spikeTimes),size(maps.rMaps{1},2)) nan(length(spikeTimes),size(maps.rMaps{1},2))}; 
for j = 1:length(spikeTimes) %  1:20 % 
    
    for itDir = 1:2
        
        % Get the map %
        map        = maps.rMaps{itDir+1}( j, : );  % Linearised rate maps are stored in an (nPosBin, nCell) array, for each direction. 
        % for CCW runs, we flip everything so directions always go left-right
        if itDir == 2
            map    = fliplr(map);
            linpos = abs(posData.linPos(:,itDir+1) - maxPixRMap); % flip positions as well
        else
            linpos = posData.linPos(:,itDir+1);
        end
       
        if ~isempty(regexp(lower(prms.trackType),'sq','once'))
            % shift map so minimum is always the last bin - this ensures
            % that we won't need to worry about field wrap around issues
            % for square track data
            [~,minInd] = min(map,[],2);
            nShift     = length(map) - minInd; %
            map        = circshift(map, nShift );
        else
            nShift     = 0;
        end
        
        % Don't run if below SI threshold or if peak FR too low
        [pkFR,pkInd]   = nanmax( map );
        if bits_per_spike(j, itDir) < prms.minSIForMapInclusion ||  pkFR < prms.minPeakRate
            continue;
        end
        
        % different modi operandi for field detection
        switch prms.fieldDetectMode
            case 'MaxArea>Thr' 
                % largest field > some rate threshold
                allFieldMask     = map >= nanmax(map(:)) .* prms.rateThrForField;
                
                fieldProps       = regionprops( allFieldMask, map, {'Area','MaxIntensity','MeanIntensity','PixelList'} );
                fieldProps       = fieldProps(  [fieldProps(:).Area] >= prms.sizeThrForField   );
                
                if isempty(fieldProps);   continue;   end % fail gracefully
                
                [~, maxFieldInd] = max( [fieldProps(:).Area] );   % Get the largest field: in future, could also try looking for highest peak or mean rates.
                fieldNumInd      = fieldProps(maxFieldInd).PixelList(:,1);
                
            case 'Pk2LocalMin'
                % Define a field on the basis of local minima flanking the largest peak %
                % 1 - find the local minima. %
                mapPad        = [1000, map, 1000];    % .. local minima on the padded map, so we will find peak flanking minima even if last or first bin.
                lclMinInd     = islocalmin( mapPad );
                lclMinInd     = lclMinInd( 2:end-1 );
                
                % 2 - remove those whose firing rate, relative to the peak, is
                % too high (unless it is an end bin, these always stay.) %
                tooHighLMInd            = [false, map(2:end-1) > (prms.localMinThrForField * pkFR), false];
                lclMinInd(tooHighLMInd) = false;
                
                % 3 - Now create a logical index of bins between those LMs that
                % flank the peak. %
                lclMinInd     = find( lclMinInd );     % Convert to numerical index.
                fStartNumInd  = find( lclMinInd < pkInd, 1, 'last' );
                fEndNumInd    = find( lclMinInd > pkInd, 1, 'first' );
                
                fieldNumInd   = lclMinInd(fStartNumInd) : lclMinInd(fEndNumInd);
                
            case 'allFields'
                % use all subfields that satisfy size and rate threshold
                mapPad        = [1000, map, 1000];
                lclMinInd     = islocalmin( mapPad );
                lclMinInd     = find( lclMinInd(2:end-1) );
                lclMinInd     = [lclMinInd(1:end-1)' lclMinInd(2:end)'];
                lclMinInd( diff(lclMinInd,[],2) < prms.sizeThrForField , : ) = []; % size thresh (too small)
                
                fieldNumInd = cell(size(lclMinInd,1),1);
                for i = 1:size(lclMinInd,1)
                    if max(map(lclMinInd(i,1):lclMinInd(i,2))) >= prms.localMinThrForField * pkFR % rate thresh
                        fieldNumInd{i} = lclMinInd(i,1):lclMinInd(i,2); % rate = high enough
                    end
                end
                fieldNumInd = fieldNumInd(~cellfun('isempty',fieldNumInd));
            
            % new case to try and deal with pup fields better using peak and trough prominence
            case 'pkProminence'
                % use all subfields with distinct and sufficicently prominent peaks

                % 1 - find the local minima and maxima. %
                mapPad        = [map(end), map, map(1)];    % .. local minima on the padded map, so we will find peak flanking minima even if last or first bin.
                [lclMinInd, lclMinProm]     = islocalmin( mapPad ); %find local minima and prominence
                [lclMaxInd, lclMaxProm]     = islocalmax( mapPad); % find local maxima and prominence

                lclMinInd     = lclMinInd( 2:end-1 ); % cut padding values
                lclMinProm    = lclMinProm( 2:end-1 );
                lclMaxInd     = lclMaxInd (2:end-1 );
                lclMaxProm    = lclMaxProm( 2:end-1 );

                % 2 - remove maxima below threshold
                for i = 1:length(map)
                    if lclMaxProm(i) < prms.maxPromThreshold % prominence
                        lclMaxProm(i) = 0;
                        lclMaxInd(i) = 0;
                    end
                    if map(i) < max(map) * prms.rateThrForField % rate
                        lclMaxProm(i) = 0;
                        lclMaxInd(i) = 0;
                    end
                end

                % 3 - threshold minima 
                for i = 1:length(map)
                    if lclMinProm(i) < prms.minPromThreshold % prominence
                        lclMinProm(i) = 0;
                        lclMinInd(i) = 0;
                    end
%                     if map(i) > max(map) * prms.rateThrForField % rate
%                         lclMinProm(i) = 0;
%                         lclMinInd(i) = 0;
%                     end
                end

                % 4 - list included maxima indices
                maxInd = find(lclMaxInd);

                % 5 - find local minima either side of the included maxima
                for j = 1:length(maxInd)
                    min2 = 0;
                    min1 = 0;
                    for i = 1:length(lclMinInd)
                        if lclMinInd(i) > 0
                            min1 = min2;
                            min2 = i;
                            if min2 > maxInd(j)
                                break
                            end
                        end
                    end
    
                    fieldNumInd{j,1} = min1 : min2;
                end
                fieldNumInd = fieldNumInd(~cellfun('isempty',fieldNumInd));
            otherwise
        end

       
        % 4 - wrap pos data (same amount as map was shifted)
        wrapPixRaw                          = nShift * prms.binSize_linMaps; % 
        wrappedPos                          = linpos + wrapPixRaw;
        wrappedPos(wrappedPos > maxPixRMap) = wrappedPos(wrappedPos > maxPixRMap) - maxPixRMap; % wrap around
        % .. also bin it 
        wrappedPosBin                       = ceil( wrappedPos ./ prms.binSize_linMaps );
        wrappedPosBin( wrappedPosBin==0 )   = 1;

        
        % 5 - now we want make sure we only use those samples that satisfy
        % these criteria :
        % a) rat is in field, 
        % b) is a run of correct direction, 
        % c) is moving at a sufficient speed

%         if ~iscell(fieldNumInd)
%             dwellInFieldInd     = ismember( wrappedPosBin, fieldNumInd );
%             speedInd            = posData.speed >= prms.speedThrForPhasePrec;
%             runDirInd           = posData.linDir == itDir;
%             validDwellSampsInd  = dwellInFieldInd  & speedInd & runDirInd;    
%         else

    % in case all sub fields are used...
        dwellInFieldInd          = false(size(posData.linPos,1),length(fieldNumInd));
        clear validSpikesInd; clear X; clear P;
        for i = 1:length(fieldNumInd) %loop though every subfield

            dwellInFieldInd(:,i) = ismember( wrappedPosBin, fieldNumInd{i} ); % accumulate index for all sub fields
        
            speedInd                 = posData.speed >= prms.speedThrForPhasePrec;
            runDirInd                = posData.linDir == itDir;
            validDwellSampsInd(:,i)       = any(dwellInFieldInd,2)  & speedInd & runDirInd; 
             
    
            % make final index... 
            spkToDwellSampInd = ceil( spikeTimes{j} .* prms.posFs );
            validSpikesInd(:,i)    = ismember( spkToDwellSampInd , find(validDwellSampsInd(:,i)) ); % 'valid' spikes are chosen on the basis position and other behavioural determinants. 
    
            % ... grab the data

            X{i}                 = wrappedPos( spkToDwellSampInd(validSpikesInd(:,i))); %#ok<*AGROW> 
            P{i}                = spikePhase{j}(validSpikesInd(:,i));

            % circular linear regression on each subfield
            [m(j,itDir,i), phi0(j,itDir,i), rho_CircReg(j,itDir,i), p_CircReg(j,itDir,i)] = circ_regress(P{i}, X{i}, prms.mNCyclForRegression, prms.minMeth);  % use Kempter et al. 2016
    %         [slope, intercept]=circRegress(X,P)
            [rho(j,itDir,i), pVal(j,itDir,i)]                                         = circ_corrcl( P{i}(~isnan(P{i})), X{i}(~isnan(P{i})) ); % not sure I understand what's the exact relationship to the rho from Kempter
    %             end
        end

        %concat all X and P fields into one double array
        Xall = [];
        Pall = [];
        for i = 1: size(X,1)
            Xall = cat(2, Xall,  cell2mat(X(i)));
            Pall = cat(2, Pall, cell2mat(P(i)));
        end

        % check theta mod
        validSpikesInd2   = ismember( spkToDwellSampInd, find(speedInd & runDirInd) ); % 'valid' spikes for whole track to estimate theta mod
        P2                = spikePhase{j}(validSpikesInd2);
        pTheta(j,itDir)   = circ_rtest(P2);
        % some more useful stuff
        nSpikesDir(j,itDir)                                                   = sum(sum(validSpikesInd));
        mapsOut{itDir}(j,:)                                                   = map;  
        allPhasePos{itDir}(length(allPhasePos{itDir})+1:length(allPhasePos{itDir})+length(Xall),:) = [Xall, Pall, j*ones(length(Xall),1) ]; % accumulate all data

        % now make a density map (phase by pos)
        binPhase                = ceil(Pall / binSzPhase);
        binPhase(binPhase == 0) = 1;
        binPos                  = ceil(((Xall - nanmin(Xall)+eps) / (nanmax(Xall)-nanmin(Xall)+eps) * prms.nPosBinsDensMap));
        densMap                 = accumarray([binPhase(~isnan(binPhase)) binPos(~isnan(binPhase))],1,[floor(360/prms.binSzPhaseDeg) prms.nPosBinsDensMap]);
        % makes sense to normalise map by spike count to even out rate differences across cells (see Cei et al. NatNeuro)
        densMapOut(:,:,j,itDir) = densMap ./ sum(densMap(:));
        
    end   
end

% output
Res.allPhaseData = allPhasePos;
Res.globalPhase0 = globalPhase0;
Res.nSpikes      = nSpikesDir;
Res.SI           = bits_per_spike;
Res.rhoCircReg   = rho_CircReg;
Res.pCircReg     = p_CircReg;
Res.circReg      = cat(4, m, phi0); %concat in 4th dim - change from original code
Res.rho          = rho;
Res.pVal         = pVal;
Res.pTheta       = pTheta;
Res.densMapRaw   = densMapOut;
Res.mapsShift    = mapsOut;
Res.rhoCircReg   = rho_CircReg;

% make table
Res                     = struct2table(Res,'AsArray',1);
Res.Properties.UserData = prms;


end


%         if strcmp( prms.fieldDetectMode, 'MaxArea>Thr' )
%             allFieldMask   = map >= nanmax(map(:)) .* prms.rateThrForField;
%             % in case any field wraps around end of track we shift map so points at end/beginning are not part of field
%             if allFieldMask(1) && allFieldMask(end)
%                 temp = strfind(allFieldMask,[0 0]);
%                 nShift = temp(1);
%                 map = circshift(map, -nShift(1));
%                 allFieldMask = circshift(allFieldMask, -nShift(1));
%             else
%                 nShift = 0;
%             end
% 
%             fieldProps     = regionprops( allFieldMask, map, {'Area','MaxIntensity','MeanIntensity','PixelList'} );
%             fieldProps     = fieldProps(  [fieldProps(:).Area] >= prms.sizeThrForField   );
% %             if isempty(fieldProps);   continue;   end
%             [~, maxFieldInd] = max( [fieldProps(:).Area] );   % Get the largest field: in future, could also try looking for highest peak or mean rates.
%             fieldNumInd = fieldProps(maxFieldInd).PixelList(:,1);
%             
%         elseif strcmp( prms.fieldDetectMode, 'Pk2LocalMin' )
%             % Define a field on the basis of local minima flanking the largest peak.
%             % 1 - find the local minima. %
% %             mapPad        = [1000, map, 1000];           % .. local minima on the padded ('wrapped') map, so we will find peak flanking minima even if they wrap the split point.
%             lclMinInd     = islocalmin( map );
% %             lclMinInd     = lclMinInd( 2:end-1 );
%             
%             % 2 - remove those whose firing rate, relative to the peak, is
%             % too high (unless it is an end bin, these always stay.) %
%             tooHighLMInd            = map > prms.localMinThrForField * pkFR;
%             lclMinInd(tooHighLMInd) = false;
%             
%             % 3 - Now create a logical index of bins between those LMs that
%             % flank the peak. %
%             lclMinInd     = find( lclMinInd );     % Convert to numerical index.
%             fStartNumInd  = find( lclMinInd < pkInd, 1, 'last' );
%             fEndNumInd    = find( lclMinInd > pkInd, 1, 'first' );
%             % in case of wrap around we need to shift things
%             if isempty(fStartNumInd)
%                 nShift = lclMinInd(end)-1;
%                 map = circshift(map, -nShift );
%                 fStart = 1;
%                 fEnd = length(map) + lclMinInd(fEndNumInd) - lclMinInd(end);
%                 
%             else
%                 nShift = 0;
%                 fStart = lclMinInd(fStartNumInd);
%                 fEnd   = lclMinInd(fEndNumInd);
%             end
%             fieldNumInd = fStart : fEnd;
%             
%         elseif strcmp( prms.fieldDetectMode, 'allFields' )
%             lclMinInd     = islocalmin(map);
%             lclMinInd     = find( lclMinInd );
%             lclMinInd     = [lclMinInd(1:end-1)' lclMinInd(2:end)'];
%             lclMinInd( diff(lclMinInd,[],2) < 5, : ) = [];
%             
%             peakThrsh = false(size(lclMinInd,2),1);
%             for i = 1:size(lclMinInd,1)
%                 peakThrsh(i) = max(map(lclMinInd(i,1):lclMinInd(i,2))) >= 0.5 * pkFR;
%             end
%             lclMinInd = lclMinInd(peakThrsh,:);
%             % check for wraparound
%             idx = find(lclMinInd(1:end-1,2) - lclMinInd(2:end,1) < -1);
%             nShift = lclMinInd(idx(1), 2);
%             map = circshift(map, -nShift );
%             lclMinInd = lclMinInd - nShift;
%             lclMinInd(lclMinInd<1) = lclMinInd(lclMinInd<1) + length(map);
%             lclMinInd = sortrows(lclMinInd);
%         end


% 
% % for end we check from other end of array
%         [~,endInd]                                  = max(fliplr(cycleIndByNSamples),[],2); 
%         endInd                                      = size(cycleIndByNSamples,2) - endInd + 1;
%         t1                                          = endInd(cycleInd); % expand
%         % get the index for each spike
%         tempL                                       = cycleLength(cycleID); % length for all legit cycles
% %         interpSampleInCycle                         = ceil(tempL(cycleInd) .* (globalPhase0Ind(~isnan(tempCycleID)) - t0) ./ (t1-t0+1)); % relative sample for spike phase within cycle
% %         interpSampleInCycle(interpSampleInCycle==0) = 1;
% %         phaseInd                                    = startInd(cycleInd) + interpSampleInCycle - 1; % actual sample - CHECK -1 IS CORRECT!!!!
% %         % 3. assign normalised phase
%         allPhase                                    = nan(length(globalPhase0Ind),1); % need to remember nan's
% %         allPhase(~isnan(tempCycleID))               = mod( eegData.phase( phaseInd ) - globalPhase0, 2*pi);
