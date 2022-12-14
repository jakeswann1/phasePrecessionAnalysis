function makePlotsPhaseAnalysis(Res,plotTypeFlag,prms,optionalPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 4
    saveFlag = true;
else
    saveFlag = false;
end

if strcmpi(plotTypeFlag,'dataset')
    ageBins = [Res.age Res.age];
    figNameStr = {[char(Res.animal) '_P' num2str(Res.age)]};
elseif strcmpi(plotTypeFlag,'population')
    ageBins = prms.ageBins;
    figNameStr = cellstr(strcat('P',num2str(ageBins(:,1)),'_',num2str(ageBins(:,2))));
else
    error([plotTypeFlag ' is not valid mate. Try ''dataset'' or ''population''.']);
end

screenSz = get(0,'ScreenSize');
for i = 1:size(ageBins,1)
    
    ageInd = Res.age >= ageBins(i,1) & Res.age <= ageBins(i,2);
    
    if sum(ageInd == 0); continue; end
    
    dataSetID = cellfun(@(x,y) y*ones(length(x),1),Res.allPhaseData(ageInd,:),num2cell([1:sum(ageInd); 1:sum(ageInd)])','uni',0);
    dataSetID = arrayfun(@(x) vertcat(dataSetID{:,x}), 1:size(dataSetID, 2), 'UniformOutput', false);
    
    % how do we define significance
    if strcmp(prms.pType,'circ_reg')
        p = vertcat(Res.pCircReg{ageInd});
    else
        p = vertcat(Res.pVal{ageInd});
    end
    
    hFig = figure('units','pixels','Position',[0.1*screenSz(3) 0.1*screenSz(4) 0.6*screenSz(3) 0.8*screenSz(4)] );
    figSz = get(hFig,'position');
    plotSz = [200 200];
    offSet = [0 0];
    
    for j = 1:2
        ax = axes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.7*figSz(4) plotSz]);
        % make population scatter
        
        phaseData4plot = [vertcat(Res.allPhaseData{ageInd,j}) dataSetID{j}];
        
        % normalise positions first
        [cellIDs,~,iB] = unique(phaseData4plot(:,3:4),'rows','stable');
        
        minPosFields = accumarray(iB,phaseData4plot(:,1),[], @min);
        maxPosFields = accumarray(iB,phaseData4plot(:,1),[], @max);
        normPos      = (phaseData4plot(:,1) - minPosFields(iB)) ./ (maxPosFields(iB) - minPosFields(iB));
        %             meanPosFields = accumarray(iB,Res.allPhaseData{j}(:,1),[],@mean);
        %             normPos       = Res.allPhaseData{j}(:,1) - meanPosFields(iB);
        %             maxPosFields  = accumarray(iB,normPos,[],@(x) max(abs(x)));
        %             normPos       = normPos./maxPosFields(iB);
        
        if prms.useSigOnly
            sigCells = find(p(:,j) < 0.05);
            sigInd   = ismember(cellIDs, sigCells);
            sigInd   = sigInd(iB);
        else
            sigInd   = true(size(phaseData4plot,1),1);
        end
        
        plot(ax,normPos(sigInd,1),repmat(phaseData4plot(sigInd,2),1,2) + [0 2*pi],'k.','MarkerSize',1);
        hold on
        [mAll, phi0All, rho_CircRegAll, p_CircRegAll,  phiOffAll] = circ_regress(phaseData4plot(sigInd,2), normPos(sigInd,1), Res.Properties.UserData.mNCyclForRegression,Res.Properties.UserData.minMeth);  % use Kempter et al. 2016
        plot(ax,[0 1]',repmat([0 1]' .* mAll + mod(phi0All,2*pi),1,3) + [0 2*pi 4*pi] ,'r-','linewidth',2);
        %             plot(ax,[-1 1]',repmat([-1 1]' .* mAll + mod(phi0All,2*pi),1,3) + [0 2*pi 4*pi] ,'r-','linewidth',2);
        hold off
        set(ax,'ylim',[0 4*pi],'xlim',[0 1],'xtick',[0 1],'ytick',[0:pi:4*pi],'yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
        %             set(ax,'ylim',[0 4*pi],'xlim',[-1 1],'xtick',[-1 1],'ytick',[0:pi:4*pi],'yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
        
        offSet(1) = offSet(1) + 2*plotSz(1) + 80;
    end
    offSet = [plotSz(1)+20 0];
    
    % make smooth density map for full data set
    for j = 1:2
        ax = axes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.7*figSz(4) plotSz]);
        
        allMapsRaw = cat(3,Res.densMapRaw{ageInd});
        
        
        if prms.useSigOnly
            sigInd   = p(:,j) < 0.05;
        else
            sigInd   = true(size(allMapsRaw,3),1);
        end
        densMapAll   = nanmean(allMapsRaw(:,:,sigInd,j),3);
        tempPad      = padarray(densMapAll,[ceil(prms.sigma(1)*2) 0],'circular');    % padd theta part of map as circular
        tempPad      = padarray(tempPad,[0 ceil(prms.sigma(2)*2)],'replicate');      % padd pos part of map
        densMapSM    = imgaussfilt(tempPad, prms.sigma,'FilterSize',prms.sigma*4+1); % smooth (filter size must be odd)
        densMapSM    = densMapSM( ceil(prms.sigma(1)*2)+1:end-ceil(prms.sigma(1)*2), ceil(prms.sigma(2)*2)+1:end-ceil(prms.sigma(2)*2) ); %t
        
        imagesc(ax,repmat(densMapSM,2,1),[0 nanmax(densMapSM(:))]);
        colormap(jet(20));
        set(ax,'ydir','normal','xtick',[],'ytick',[]);
        
        offSet = offSet(1) + 2*plotSz(1) + 80;
    end
    offSet = [0 0];
    
    % plot distribution of slopes and phi's
    for j = 1:2
        
        sigInd   = p(:,j) < 0.05;
        
        circRegData = vertcat(Res.circReg{ageInd});
        
        pax = polaraxes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.4*figSz(4) 0.9*plotSz]);
        polarhistogram(pax,circRegData(sigInd,j,2),linspace(15*pi/180,2*pi+15*pi/180,13),'facecolor','r');
        hold(pax,'on')
        polarhistogram(pax,circRegData(~sigInd,j,2),linspace(15*pi/180,2*pi+15*pi/180,13),'facecolor',[0.5 0.5 0.5]);
        hold(pax,'off')
        
        offSet(1) = offSet(1) + plotSz(1) + 20;
        
        ax = axes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.4*figSz(4) plotSz]);
        if sum(sigInd) > 0
            plot(ax,[0 1]',[0 1]' .* circRegData(sigInd,j,1)' ,'k-','linewidth',0.5);
        end
        hold(ax,'on');
        plot(ax,[0 1]',[0 1]' .* circRegData(~sigInd)' ,'linestyle',':','color',[0.5 0.5 0.5],'linewidth',0.5);
        plot([0 1]',[0 1]' .* nanmean(circRegData(sigInd,j,1)),'r-','linewidth',3);
        plot([0 1]',[0 1]' .* nanmean(circRegData(~sigInd,j,1)),'r--','linewidth',1.5);
        plot([0 1]',[0 0]','k--');
        hold(ax,'off');
        set(ax,'xlim', [0 1],'xtick',[],'ytick',0);
        
        
        offSet(1) = offSet(1) + plotSz(1) + 65;
    end
    offSet = [0 0];
    
    for j = 1:2
        ax = axes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.1*figSz(4) plotSz]);
        bar(ax,[sum(p(:,j) < 0.05)/sum(~isnan(p(:,j))); sum(p(:,j) >= 0.05)/sum(~isnan(p(:,j)))]);
        set(ax,'xticklabel',{'sig', 'ns'},'ylim',[0 0.8]);
        offSet(1) = offSet(1) + plotSz(1) + 20;
        
        ax = axes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.1*figSz(4) plotSz]);
        
        allPhaseData = vertcat(Res.allPhaseData{ageInd,j});
        
        bins = linspace(0,2*pi,361);
        rateHist = histcounts(allPhaseData(:,2), bins );
        rateHist = imfilter(rateHist,fspecial('gaussian',[1 180],45),'circular'); % smooth with Gaussian (sigma=45deg, so we get very smooth curve)
        plot(rateHist,'b-','linewidth',2);
        set(ax,'xlim',[0 360],'xtick',[0,180,360],'xticklabel',{'0' '\pi' '2\pi'},'ylim',[min(rateHist) max(rateHist)]);
        offSet(1) = offSet(1) + plotSz(1) + 65;
    end
    
    
%     % save & close
%     if saveFlag
%         saveFigAsPDF( [char(figNameStr{i}) '_allCellsPooled'], [fileparts(optionalPath) filesep] );
%         close(gcf);
%     end
    
    if prms.plotAllCells
        
        % plot all cells in dataset
        offSet   = [20 20];
        plotSep  = [20 20];
        plotSize = [200 100];
        
        hScroll = scanpix.plot.createScrollPlot( [0.1*screenSz(3) 0.1*screenSz(4) 3*plotSize(1)+8*plotSep(1) 0.8*screenSz(4) ] );
        
        mapsAll      = { vertcat(Res.mapsShift(ageInd,1)), vertcat(Res.mapsShift(ageInd,2)) };
        allPhaseData = { vertcat(Res.allPhaseData(ageInd,1)), vertcat(Res.allPhaseData(ageInd,2)) };
        allCircReg   = vertcat(Res.circReg(ageInd));
        
        for m = 1:length(ageInd)
            cellIt = cellIDs(cellIDs(:,2)==m,1);
            for j = cellIt' % cellID
                
                for k = 1:2 % directions %

                    map2plot = mapsAll{k}{m}(j,:);
                    ind = allPhaseData{k}{m}(:,3) == j;
                    if sum(ind) == 0; continue; end
                    
                    temp_P = allPhaseData{k}{m}(ind,2); % phase for cell
                    temp_P = repmat(temp_P,1,2) + [0 2*pi];
                    temp_X = allPhaseData{k}{m}(ind,1); % pos for cell
                    
                    % PLOT LEFT HAND PLOT (WHOLE FIELD)
                    hAx = scanpix.plot.addAxisScrollPlot( hScroll, [offSet plotSize], plotSep );
                    area(hAx,map2plot); % plot map
                    hold(hAx,'on');

                    % mark field with dashed vertical lines
                    plot(hAx,[min(Res.Xsub{1}{j}{1}./Res.Properties.UserData.binSize_linMaps) min(Res.Xsub{1}{j}{1}./Res.Properties.UserData.binSize_linMaps);max(Res.Xsub{1}{j}{1}./Res.Properties.UserData.binSize_linMaps) max(Res.Xsub{1}{j}{1}./Res.Properties.UserData.binSize_linMaps)]',...
                        [0 ceil(nanmax(map2plot));0 ceil(nanmax(map2plot))]','k--');
                    hold(hAx,'off');
                    set(hAx,'xtick',[],'ylim',[0 ceil(nanmax(map2plot))]);
                    offSet(1) = offSet(1) + plotSize(1) + plotSep(1);

                    % ########### Edited henceforth to account for multiple
                    % possible subfields ###########

                    % rearrage regressions so each cell has subfields along
                    % the first dimension and remove zero values
                    fields = size(allCircReg{m});
                    fieldsFirst = permute(allCircReg{m},[3 4 1 2]);
                    cellFields = mat2cell(fieldsFirst, [fields(3)], [fields(4)], [ones(fields(1),1)], [ones(fields(2),1)]);
                    cellFields = cellfun(@(c) c(c~=0), cellFields, 'UniformOutput', false);

                    for x = 1:fields(1)
                        for y = 1:fields(2)
                            cellFields{:,:,x,y} = transpose(reshape(cellFields{:,:,x,y}, [],2));
                        end
                    end
                    % Regression data are now structured: 4-D cell array with:
                    % dim1&2 = regression values for each subfield, values
                    % for an individual subfield in dim1
                    %dim3 = cell (j)
                    %dim4 = direction (k)

                    %PLOT RIGHT HAND SIDE PLOT
                    hAx = scanpix.plot.addAxisScrollPlot( hScroll, [offSet plotSize(1)/2 plotSize(2)], plotSep );
                    
                    % plot(hAx,temp_X./Res.Properties.UserData.binSize_linMaps, temp_P./(4*pi), 'k.'); 
                    % plot spikepos v phase in field 
                    hold(hAx,'on');

                    % plot rate vs position RHS
                    plot(hAx,map2plot./nanmax(map2plot(:)),'b-','linewidth',2); % plot field

%                     try
                        for n = 1:size(Res.Xsub{1}{j},2) % subfields

                            % calculate circular regression line values
                            % NB currently restricted to negative
                            % regression values only
                            yVals{n,j,k} = repmat([0 1]' .* allCircReg{m}(j,k,n,1) + mod(allCircReg{m}(j,k,n,2),2*pi),1,3) + [0 2*pi 4*pi];
                            %yVals = repmat([-1 1]' .* Res.circReg{1}(j,k,1) + mod(Res.circReg{1}(j,k,2),2*pi),1,3) + [0 2*pi 4*pi]; % plot circular regression line
                        
    
                            %plot spikes vs phase for the current subfield
                            plot(hAx,Res.Xsub{1}{j}{n}./Res.Properties.UserData.binSize_linMaps, Res.Psub{1}{j}{n}.*(2/(4*pi)), 'k.');
    
                            %set x boundaries per subfield
                            rangeTop = max(Res.Xsub{1}{j}{n})./Res.Properties.UserData.binSize_linMaps;
                            try
                                rangeBot = min(setdiff(Res.Xsub{1}{j}{n}, Res.Xsub{1}{j}{n-1}))./Res.Properties.UserData.binSize_linMaps;
                            catch
                                rangeBot = min(Res.Xsub{1}{j}{n})./Res.Properties.UserData.binSize_linMaps;
                            end
                            
                            %plot dashed lines for P>0.05, unbroken line for
                            %P<0.05
                            if p(j,k,n) >= 0.05
                                plot(hAx,[rangeBot rangeTop],yVals{n,j,k}/(2*pi),'r--');
                                hold(hAx, 'on');
    %                             subtitle('p = '+string(p(j,k)))

                            else
                                plot(hAx,[rangeBot rangeTop],yVals{n,j,k}/(2*pi),'g-');
                                hold(hAx,'on');
    %                             subtitle('p = '+string(p(j,k)))
                            end
                        
                        end
                    
                    % set axes
                    set(hAx,'xlim',[0.95*min(Res.Xsub{1}{j}{n}./Res.Properties.UserData.binSize_linMaps) 1.05*max(Res.Xsub{1}{j}{n}./Res.Properties.UserData.binSize_linMaps)],'ytick',[],'xtick',[],'ylim',[-0.05 1]);
                    offSet(1) = offSet(1) + plotSize(1)/2 + 3*plotSep(1);

%                     catch
% %                         PLOT LEFT HAND PLOT (WHOLE FIELD)
% %                         map2plot = mapsAll{k}{m}(j,:);
% %                         hAx = scanpix.plot.addAxisScrollPlot( hScroll, [offSet plotSize], plotSep );
% %                         area(hAx,map2plot); % plot map
% %                         hold(hAx,'on');
% 
%                         set(hAx,'xlim',[0 100],'ytick',[],'xtick',[],'ylim',[-0.2 1]);
%                         offSet(1) = offSet(1) + plotSize(1)/2 + 3*plotSep(1);
% 
%                     end
                end
                offSet(1) = 20;
                offSet(2) = offSet(2) + plotSize(2) + plotSep(2);
                hold(hAx, 'off');
                
            end
        end


        end

%         % save & close
%         if saveFlag
%             saveFigAsPDF( [char(figNameStr{i}) '_singleCells'],[fileparts(optionalPath) filesep]);
%             close(gcf);
%         end
    end
    
    
    %         case 'population'
    %
    %             screenSz = get(0,'ScreenSize');
    %             hFig = figure('units','pixels','Position',[0.1*screenSz(3) 0.1*screenSz(4) 0.6*screenSz(3) 0.8*screenSz(4)] );
    %             figSz = get(hFig,'position');
    %             plotSz = [200 200];
    %             offSet = [0 0];
    %
    %             for i = 1:2
    %
    %                 for j = 1:size(prms.ageBins,1)
    %
    %
    %
    %                     ax = axes('units','pixels','position',[0.05*figSz(3)+offSet(1) 0.7*figSz(4) plotSz]);
    %
    %
    %                 end
    %
    %
    %
    %
    %             end
    %
    %
    %
    %
    %
    %         otherwise
    %             error([plotTypeFlag ' is not valid mate. Try ''dataset'' or ''population''.']);
    %
    %     end
end


