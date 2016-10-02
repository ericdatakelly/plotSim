function [] = plotSim(simulationName)

% plotSim.m
% syntax: plotSim('rockName_simulationName')
% Program for plotting a summary of CRYSTALLIZE3D simulation results.
%
% This version uses the formatting of input and output files from
% Crystallize3D v1.4.
% The following files must be present in the subdirectory 'RockData':
%   rockName.txt
% The following files must be present in the subdirectory 'SimulationData':
%   params_simulationName.txt
%   simulationName.log
%   simulationName.int.txt
%   simulationNamerun.txt
%   simulationNamets#####xl.txt (the xl file with the largest time step #####)
% The following file must be present in the subdirectory 'PlottingParams':
%   rockNamePlottingParams.txt
% The first file is the CT data (be sure there are 5 header lines if the int
% file was created before xllize1.1), the second is the Crystallize parameters file,
% the next four are output by Crystallize, and the last is a user-created file
% of rock data and plotting parameters.  If no rock plotting params file is
% found, a generic one will be used instead.
%
% The file startup.m is used to change the directory to the plotSim
% directory.  Be sure to update it as needed.  It is typically installed
% here: \toolbox\local

%IMPROVEMENTS:
% Move font, line, margin properties to a txt file - make one file that is
%   used for all samples.
% Use a text file for other options like the plotSim working directory, the
%   plotting params directory, the data directory, etc.
% How can I make the file naming convention better?
%   Prompt the user to select a rock name(or no rock).
%   Default rock name is the last name.
%   If no rock name, skip the plotting functions for the rock data.  Use a
%       GUI to show a list of models available (with the newest model at
%       the top).
%   Use 'fullfile' to build a path and filename more efficiently.
%   'readfile' can search for parameters in a file.  Use this to search the
%       params file for the parameters instead of designating exact line
%       locations.  Updates to Crystallize3D tend to force me to recode
%       plotSim to fit the new params file.  This means that plotSim is not
%       backwards compatible with older Crystallize output.  On that topic,
%       consider what to do with parameters that are obsolete (maybe just
%       allow plotSim to read and store them even if we don't use them
%       anymore).
% Use separate m files to create each figure (3D plot, layers plot, and
% main plot.
% Many fid assignments should be closed if the file does not exist (close
% them before 'return')


% Set up paths and constants
fullPath = mfilename('fullpath');
plotSimIndex = strfind(fullPath,'plotSim');
plotSimPath = fullPath(1:(plotSimIndex(1)-1));
rockDataPath = strcat(plotSimPath,'RockData\');
simulationDataPath = strcat(plotSimPath,'SimulationResults\');
paramsPath = strcat(plotSimPath,'PlottingParams\');
R = 0.0083144; % Gas constant in kJ/mol/K

% Set some of the font, line, and color characteristics
fontSizeForStack = 10;
lineWidthForStack = 1;
fontSizeForCSD = 10;
lineWidthForCSD = 2;
barWidthForCSD = 0.7;
lengthOfTicks = [0.02 0.02]; % [2D 3D] fraction of axis length
fontSizeRxtDistribution = 10;
sizeOfText = 8; % Font size for text lists
baseColor = 0.2; % This sets the rxt dist plot background value (grayscale between 0 and 1)

% Set the positions of the figures and other positions within figures
% Main figure position (Sim2 is the identifier for this figure)
screenSize = get(0,'ScreenSize');
figPosWidthSim2 = (3/8)*screenSize(3);
figPosheightSim2 = (11/8.5)*figPosWidthSim2;
figPosBottomSim2 = screenSize(4)-figPosheightSim2-85;
figPosLeftSim2 = -40;
figPositionSim2 = ...
    [figPosLeftSim2 figPosBottomSim2 figPosWidthSim2 figPosheightSim2];

% Other positions within figures [left bottom width height]
rxtPositionSim2 = [0.100 0.600 0.350 0.320]; % Reactant distribution
csdPositionSim2 = [0.600 0.750 0.350 0.200]; % CSD
plotStackPositionSim2 = [0.700 0.070]; % Nucl rates, DrG, heating path, Al in system
plotStackDimensionsSim2 = [0.250 0.180]; % Width and height of largest plot
textPositionSim2 = [0.050 0.050]; % [left margin, bottom margin]
textIncrementSim2 = 0.018; % Spacing between lines of text
textGapSim2 = 0.010; % Gap between blocks of text
numberOfLinesSim2 = 28; % Lines of text to print
xDescPositionSim2 = 0.050; % Left margin of text descriptions
xSimPositionSim2 = 0.380; % Left margin of simulation values
xRockPositionSim2 = 0.500; % Left margin of rock values
titleRockPositionSim2 = [0.050 0.970];
titleSimPositionSim2 = [0.050 0.940];

% Read the rock-specific plotting data
if strfind(simulationName,'_');
    rockNameIndex = strfind(simulationName,'_'); % Find the rock name
    rockName = simulationName(1:rockNameIndex(1)-1); % Store the rock name
    % Check for rock files, use generic files otherwise
    % These cases are a temporary fix.
    % Instead of these cases, scan the plotting parameters directory for a
    % file with the rock name. If it exists, rockPresent == 1, etc.
    switch rockName
        case 'PM1'
            rockPresent = 1;
        case 'PM2'
            rockPresent = 1;
        case 'PM3'
            rockPresent = 1;
        case 'PM4'
            rockPresent = 1;
        case '711A'
            rockPresent = 1;
        case '160A'
            rockPresent = 1;
        case '191A'
            rockPresent = 1;
        case 'AG4'
            rockPresent = 1;
        case 'MD'
            rockPresent = 1;
        case 'WR1bt'
            rockPresent = 1;
        case 'WR1tp'
            rockPresent = 1;
        case 'WR3m'
            rockPresent = 1;
        case 'Jen'
            rockPresent = 1;
        case 'HE1'
            rockPresent = 1;
        otherwise
            rockName = 'NoName';
            rockPresent = 0;
    end
else
    rockName = 'NoName';
    rockPresent = 0;
end


% Read the plotting parameters file
filenamePlottingParamsRock = strcat(paramsPath,rockName,'PlottingParams.txt');
fid = fopen(filenamePlottingParamsRock);
if fid == -1
    disp('---------------------------------------------');
    disp('   Error in plotSim.m');
    disp('   Cannot find plotting parameters file');
    disp(['   "',strcat(rockName,'PlottingParams.txt'),'" is not in the directory']);
    disp('---------------------------------------------');
    return
else
    fclose(fid);
end
[rockFileName,grtZoningFileName,radiusVsCoreConcFileName,crystalNumberDensityRock,...
    timeRangeString,temperatureEquilibriumRock,temperature95Rock,...
    rockMode,rockVol,binMin,binCenter,binMax,csdXLimMin,csdXLimMax,csdXTickMin,...
    csdXTickIncrement,csdXTickMax,csdYLimMin,csdYLimMax,plotLayerCSDs,...
    stackXLimMin,stackXLimMax,nucleationRateYAxis,nucleationRateYLimMin,...
    nucleationRateYTickInc,nucleationRateYLimMax,smoothingSwitch,...
    smoothingWindowNucl,heatingPathYLimMin,heatingPathYLimInc,...
    heatingPathYLimMax,plotCrystalsIn3D,affinityYAxis,...
    affinityYLimMin,affinityYLimInc,affinityYLimMax,plotZoning,...
    coreConcXMin,coreConcXMax,coreConcBinMin,coreConcBinCenter,...
    coreConcBinMax,smoothingWindowZoning...
    ] = readPlottingParamsRock(filenamePlottingParamsRock);

% Read the crystal radius data from the natural sample int file
if rockPresent == 1
    filenameIntRock = strcat(rockDataPath,rockFileName);
    fid = fopen(filenameIntRock);
    if fid == -1
        disp('---------------------------------------------');
        disp('   Error in plotSim.m');
        disp('   Cannot find natural sample int file');
        disp(['   "',rockFileName,'" is not in the directory']);
        disp('---------------------------------------------');
        return
    else
        fclose(fid);
    end
    [rockValuesInt] = readIntFile(filenameIntRock);
end

% Read the parameters file
filenameParams = strcat(simulationDataPath,'params_',simulationName,'.txt');
fid = fopen(filenameParams);
if fid == -1
    disp('---------------------------------------------');
    disp('   Error in plotSim.m');
    disp('   Cannot find simulation parameters file');
    disp(['   "',strcat('params_',simulationName,'.txt'),'" is not in the directory']);
    disp('---------------------------------------------');
    return
else
    fclose(fid);
end
[reportInterval,porosity,tortuosity,Q,Dinf,Teq,CiAl,k2,k1,max_ppb,max_mode,...
    voxelEdgeLength,timeStep,numVoxelsX,numVoxelsY,numVoxelsZ,...
    defaultRxtAmount,numStructures,hetType,structures,fv,n_GtoC...
    ] = readParams(filenameParams);

% Read the output files from the CRYSTALLIZE3D simulation
filenameSimulationInt = strcat(simulationDataPath,simulationName,'.int.txt');
fid = fopen(filenameSimulationInt);
if fid == -1
    disp('---------------------------------------------');
    disp('   Error in plotSim.m');
    disp('   Cannot find simulation int file');
    disp(['   "',strcat(simulationName,'.int.txt'),'" is not in the directory']);
    disp('---------------------------------------------');
    return
else
    fclose(fid);
end
[simulationValuesInt] = readIntFile(filenameSimulationInt);
% Check that there are crystals in the simulation.  When a model
% extends beyond the heating path, it will quit, and in some cases will
% quit before the first nucleation.
if size(simulationValuesInt) == 0
    disp('---------------------------------------------');
    disp('   Error in plotSim.m');
    disp('   Simulation int file contains no crystals');
    disp('---------------------------------------------');
    return
end

filenameSimulationRun = strcat(simulationDataPath,simulationName,'run.txt');
fid = fopen(filenameSimulationRun);
if fid == -1
    disp('---------------------------------------------');
    disp('   Error in plotSim.m');
    disp('   Cannot find simulation run file');
    disp(['   "',strcat(simulationName,'run.txt'),'" is not in the directory']);
    disp('---------------------------------------------');
    return
else
    fclose(fid);
end
[simulationValuesRun] = readSimRunFile(filenameSimulationRun);

% Find the final time step to determine the name of the Xl file
% The time step is a five digit string (it might start with zeros)
simulationDuration = num2str(simulationValuesRun...
    (length(simulationValuesRun(:,1)),1),'%05d');

filenameSimulationXl = strcat(simulationDataPath,simulationName,'ts',...
    simulationDuration,'xl.txt');
fid = fopen(filenameSimulationXl);
if fid == -1
    disp('----------------------------------------');
    disp('   Cannot find xl file');
    disp('   Simulation mode not calculated');
    disp('   Zoning vs Radius relationship not calculated');
    disp('----------------------------------------');
    simulationValuesXl = zeros(1,8);
    fclose(fid);
else
    fclose(fid);
    [simulationValuesXl] = readSimXlFile(filenameSimulationXl);
end

if plotZoning
    filenamSimulationValuesRadiiVsTime = strcat(simulationDataPath,simulationName,...
        '_radii-vs-time.txt');
    fid = fopen(filenamSimulationValuesRadiiVsTime);
    
    % Check that the radii-vs-time file is present
    if fid == -1
        disp('.....................................................');
        disp('   Cannot find simulation radii-vs-time file');
        disp(['   "',strcat(simulationName,'_radii-vs-time.txt'),...
            '" is not in the directory']);
        disp(' ')
        disp('   Dependent functions will not be run.')
        disp('   (Zoning-radius relationships will not be plotted)')
        disp('.....................................................');
        
        % If the file is not present, do not try to use radii-vs-time
        % values anywhere else
        plotZoning = 0;
    else
        % If the file is present, continue extracting values and checking
        % for other related files.
        fclose(fid);
        [simulationValuesRadiiVsTime] = ...
            readSimRadiiVsTimeFile(filenamSimulationValuesRadiiVsTime);
        
        % Get the radii vs time values for the first crystal
        simulationFirstRadiiVsTime = (simulationValuesRadiiVsTime(1:2,:))';
                
        filenameRadiusVsCoreConcFileName = strcat(rockDataPath,radiusVsCoreConcFileName);
        fid = fopen(filenameRadiusVsCoreConcFileName);
        if fid == -1
            disp('---------------------------------------------');
            disp('   Error in plotSim.m');
            disp('   Cannot find core conc vs radius file');
            disp(['   "',radiusVsCoreConcFileName,'" is not in the directory']);
            disp('---------------------------------------------');
            return
        else
            fclose(fid);
            [rockValuesRadiusVsCoreConc] = ...
                readRockZoningFile(filenameRadiusVsCoreConcFileName);
        end
        
        filenameGrtZoningProfileRock = strcat(rockDataPath,grtZoningFileName);
        fid = fopen(filenameGrtZoningProfileRock);
        if fid == -1
            disp('---------------------------------------------');
            disp('   Error in plotSim.m');
            disp('   Cannot find garnet zoning file');
            disp(['   "',grtZoningFileName,'" is not in the directory']);
            disp('---------------------------------------------');
            return
        else
            fclose(fid);
            [rockValuesGrtZoning] = readRockZoningFile(filenameGrtZoningProfileRock);
        end
    end
end


% Create the vectors for the plot axes limits
bins = (binMin:binCenter:binMax);                                           % CSD bins [min:center:max]
csdXLim = [csdXLimMin csdXLimMax];                                          % CSD X axis [min max]
csdXTick = (csdXTickMin:csdXTickIncrement:csdXTickMax);                     % CSD X axis tick locations [min:increment:max]
csdYLim = [csdYLimMin csdYLimMax];                                          % CSD Y axis [min max]
stackXLim = [stackXLimMin stackXLimMax];                                    % Nuclei and Nucleation Rate X axis limits [min max]
nucleationRateYLim = [nucleationRateYLimMin nucleationRateYLimMax];         % Nucleation rate Y axis limits [min max]
nucleationRateYTick = (nucleationRateYLimMin:nucleationRateYTickInc:nucleationRateYLimMax); % Nucleation rate Y axis tick marks [min:increment:max]
heatingPathYlim = [heatingPathYLimMin heatingPathYLimMax];                  % Heating Path Y axis limits [min max]
heatingPathYTick = (heatingPathYLimMin:heatingPathYLimInc:heatingPathYLimMax); % Heating Path Y axis tick marks [min:incrememnt:max]
affinityYLim = [affinityYLimMin affinityYLimMax];                           % Reaction affinity Y axis limits [min max]
affinityYTick = (affinityYLimMin:affinityYLimInc:affinityYLimMax);          % Reaction affinity Y axis limits [min:increment;max]
percentAlYLim = [0 100];                                                    % Percent Al in System Y axis limits [min max]
percentAlYTick = (25:25:100);                                               % Percent Al in System Y axis tick marks [min:increment:max]
zoningVsRadiusXLim = [coreConcXMin coreConcXMax];                           % Garnet chemical zoning vs radius X axis limits
binsCoreConc = (coreConcBinMin:coreConcBinCenter:coreConcBinMax);           % Garnet chemical zoning vs radius bins

% Remove first and last tick labels for some axes so the plots look nice
heatingPathYTick = heatingPathYTick(2:length(heatingPathYTick) - 1);
percentAlYTick = percentAlYTick(1:length(percentAlYTick)-1);

% Calculate Simulation Volume
[simulationVolume] = calcVolume(voxelEdgeLength,numVoxelsX,numVoxelsY,...
    numVoxelsZ);

% Calculate Crystal Size Distribution values
[simulationElements] = calcCSD(simulationValuesInt(:,5),bins);
if rockPresent == 1
    [rockElements] = calcCSD(rockValuesInt(:,5),bins);
end

% Calculate Number of Crystals at each time step
[timeMy,crystalsPerSecond] = ...
    calcCrystalsPerTime(simulationVolume,simulationValuesRun,...
    smoothingSwitch,smoothingWindowNucl);

% Calculate Al in reactants, products, fluid
[reactantAl,productAl,fluidAl] = calcAl(simulationValuesRun);

% Calculate Nucleation rate through time
[temperature,nucleationRateMax] = calcNuclRate(simulationValuesRun,CiAl,...
    k1,k2,n_GtoC,Teq);

% Calculate effective diffusivity at the characteristic temperature, D(Tc).
% D(Tc) is the diffusivity that will produce the same integrated
% diffusivity for the duration of the event if D(Tc) is fixed for the
% duration of the event.
[timeD,tempD,D] = calcD(Dinf,porosity,tortuosity,Q,R,simulationValuesRun);

% Calculate Min, Max, and Mean radius of cystals in rock and simulation
if rockPresent == 1
    rockRadiiMin = min(rockValuesInt(:,5));
    rockRadiiMean = mean(rockValuesInt(:,5));
    rockRadiiMax = max(rockValuesInt(:,5));
end
simulationRadiiMin = min(simulationValuesInt(:,5));
simulationRadiiMean = mean(simulationValuesInt(:,5));
simulationRadiiMax = max(simulationValuesInt(:,5));

% Calculate Mode of crystals in simulation
[simulationPorphMode] = calcMode(simulationVolume,simulationValuesXl(:,6));

% Calculate Mode of crystals specified by user (convert from vol frac of
% CAP to vol% of ppb)
max_mode_ppb = 100 * max_mode * fv;

% Calculate crystal number density
[crystalNumberDensitySim] = calcCrystalNumDensity(simulationValuesRun,...
    simulationVolume);

% Calculate the index of the first and the last nucleation (T95)
[firstCrystalIndex,lastCrystalIndex95] = calcFirstLastNucl(simulationValuesRun);

% Calculate Time range of cystallization
timeFirstCrystal = timeMy(firstCrystalIndex);  % To nearest reporting interval
timeCrystal95 = timeMy(lastCrystalIndex95); % To nearest reporting interval
timeRange = timeCrystal95 - timeFirstCrystal;

% Calculate Temperature of first and last crystals, and thermal overstepping
temperatureStart = simulationValuesInt(1,7); % This is not dependent on reporting interval (it is from the int file)
temperature95 = temperature(lastCrystalIndex95); % To nearest reporting interval
delaT = temperatureStart - Teq; % This is not dependent on reporting interval (it is from the int file)

% Calculate the first, last, and mean reaction affinity values
affinityNucl = simulationValuesXl(:,8);
affinityMean = simulationValuesRun(:,11);
affinityFirstNuclIndex = find(affinityNucl,1);
affinityNuclFirst = affinityNucl(affinityFirstNuclIndex);
affinityNuclFirstTime = (simulationValuesXl(affinityFirstNuclIndex,1)*10^(-6)); % For plotting first nucleation
affinityNuclMax = max(affinityNucl);
affinityMeanMax = max(affinityMean);

% Calculate values for plotting the Zoning vs Radius relationship
if plotZoning
    [simulationRadiusVsCoreConc,rockValuesGrtZoningNorm,...
        rockZoningFromRadiusSmoothed,simulationTimeAndZoning,...
        simulationTimeAndZoningSplineCoef] = ...
        calcZoningVsRadius(simulationValuesXl,rockValuesGrtZoning,...
        simulationFirstRadiiVsTime,smoothingWindowZoning);
    
    % Calculate differences in radius range and crystal volume at specified
    % core concentrations (bins).
    
    % Make bins to organize the results into discrete points.  Ensure that
    % the range of bins encompasses both the rock and the simulation values
    % so that they can be compared directly.
    maxConc = max(rockValuesRadiusVsCoreConc(:,2));
    if max(simulationRadiusVsCoreConc(:,2)) > maxConc
        maxConc = max(simulationRadiusVsCoreConc(:,2));
    end
    
    % Calculate the values from the rock and the simulation
    [binnedRadiusRangeAndVolRock] = ...
        calcRadRangeVol(rockValuesRadiusVsCoreConc,binsCoreConc,rockVol);
    [binnedRadiusRangeAndVolSim] = ...
        calcRadRangeVol(simulationRadiusVsCoreConc,binsCoreConc,...
        simulationVolume);
    
    % Write the values too...
    fid = fopen(strcat(simulationName,'_coreConcBins.txt'),'w');
    fprintf(fid,'Conc bins\n');
    for i = 1:length(binsCoreConc)
        fprintf(fid,'%g\n',binsCoreConc(i));
    end
    fclose(fid);
    fid = fopen(strcat(simulationName,'_rockRadVsCoreConc.txt'),'w');
    fprintf(fid,'rock r\trock conc\n');
    for i = 1:length(rockValuesRadiusVsCoreConc(:,1))
        fprintf(fid,'%g\t%g\n',...
            rockValuesRadiusVsCoreConc(i,1),...
            rockValuesRadiusVsCoreConc(i,2)...
            );
    end
    fclose(fid);
    fid = fopen(strcat(simulationName,'_simRadVsCoreConc.txt'),'w');
    fprintf(fid,'sim r\tsim conc\n');
    for i = 1:length(simulationRadiusVsCoreConc(:,1))
        fprintf(fid,'%g\t%g\n',...
            simulationRadiusVsCoreConc(i,1),...
            simulationRadiusVsCoreConc(i,2)...
            );
    end
    fclose(fid);

    
    % Calculate the differences between the rock and sim values and then
    % remove bins with no differences in radius or no volume.
    
    radiusRangeRock = binnedRadiusRangeAndVolRock(:,3) - ...
        binnedRadiusRangeAndVolRock(:,2);
    radiusRangeSim = binnedRadiusRangeAndVolSim(:,3) - ...
        binnedRadiusRangeAndVolSim(:,2);
    
    culledBins = [];
    culledRadiusRangeRock = [];
    culledRadiusRangeSim = [];
    ratioOfVol = [];
    for i = 1:length(binnedRadiusRangeAndVolRock(:,1))
        if (radiusRangeRock(i) > 0) && (radiusRangeSim(i) > 0)
            culledBins = [culledBins; binnedRadiusRangeAndVolRock(i,1)];
            culledRadiusRangeRock = [culledRadiusRangeRock; radiusRangeRock(i)];
            culledRadiusRangeSim = [culledRadiusRangeSim; radiusRangeSim(i)];
            ratioOfVol = [ratioOfVol; binnedRadiusRangeAndVolSim(i,4) /...
                binnedRadiusRangeAndVolRock(i,4)];
        end
    end
    ratioOfRadiusRange = culledRadiusRangeSim ./ culledRadiusRangeRock;
end

% Determine which crystals are in which layers for plotting CSDs from
% individual layers
if plotLayerCSDs
    if (hetType && (numStructures > 0))
        [layerBounds,crystalsInLayers] = calcCrystalsInLayers(numVoxelsX,...
            voxelEdgeLength,simulationValuesInt,structures);
    end
end


% Plot the CSDs from each layer
% There is not much room on a figure and most models don't have more
% than eight layers so this plots a maximum of eight layers.  Most models
% are also symmetrical so eight layers should contain the representative CSDs

if plotLayerCSDs
    if hetType && (numStructures > 0)
        figure('position',figPositionSim2,...
            'PaperPosition',[0.25,0.25,8,10.5],...
            'FileName',strcat(simulationName,'_Plot2'))
        % [left margin, bottom margin, width, height]
        % Adding the filename here allows the model name to be the default in
        % the Save As dialog.
        %     set(gcf,'PaperPosition',[0.25,0.25,8,10.5]);
        
        numPlots = length(crystalsInLayers(1,:));
        if numPlots > 8
            numPlots = 8;
        end
        for i=1:numPlots
            indices = find(crystalsInLayers(:,i));
            h_layerCSD_axes = subplot(4,2,i);  % (rows columns plot-number(left-to-right))
            if find(indices) % if the layer contains crystals
                numCrystalsInLayers = hist(crystalsInLayers(indices,i),bins); % can't use calcCSD because I need to normalize by the total number of crystals
                numCrystalsInLayers = 100*(numCrystalsInLayers./length(simulationValuesInt(:,5)));
                if rockPresent == 1  % If the natural sample is present, plot the CSD blue line for comparison with the simulation
                    %h_layerCSD = plot(bins,rockElements,'b','LineWidth',3);
                    plot(bins,rockElements,'b','LineWidth',3);
                    hold on
                end
                bar(bins,numCrystalsInLayers,0.7,'y')
                hold off
                set(h_layerCSD_axes,...
                    'FontSize',10,...
                    'TickDir','out',...
                    'XTick',csdXTick,...
                    'XLim',csdXLim,...
                    'YLim',csdYLim);
                title(['X Voxels ',num2str(layerBounds(i,1)),'-',num2str(layerBounds(i,2))])
                if (rockPresent == 1) && (i == 1)
                    legend('Rock','Simulation','Location','Best')
                end
                if i > 6
                    xlabel('Radius (cm)')
                end
                ylabel('Crystals (%)')
            else
                bar(0,0)
                title(['X Voxels ',num2str(layerBounds(i,1)),...
                    '-',num2str(layerBounds(i,2)),' (No Crystals)'])
            end
        end
        
        % Title
        axes('position',[0 0 1 1],'visible','off') % Must establish axes in order to write text to a figure
        % Position = [left bottom width height]
        if rockPresent == 1
            text(0.150,0.970,rockName,'FontWeight','Bold','FontSize',16)
        end
        text(0.550,0.970,['Simulation: ',simulationName],'Interpreter','none','FontWeight','Bold','FontSize',12)
    else
        fprintf('\nCSDs for layers cannot be plotted because')
        fprintf('the model does not contain layers.  To remove')
        fprintf('this message, turn off the plotting option.\n')
    end
end


% Plot the remaining results in a new figure

% Make a figure on which to plot everything
% Position = [left bottom width height]
figure('Position',[figPositionSim2(1)+50 figPositionSim2(2) figPositionSim2(3) figPositionSim2(4)],...
    'PaperPosition',[0.25,0.25,8,10.5],...
    'FileName',strcat(simulationName,'_Plot1'))
%     set(gcf,'PaperPosition',[0.25,0.25,8,10.5]); %[left margin, bottom margin, width, height]

% Add Reactant Distribution
h_rxt_axes = axes;
% Define color map
colorMapGreen = (baseColor:0.005:1)';
colorMapRed = zeros(size(colorMapGreen)) + baseColor;
colorMapBlue = zeros(size(colorMapGreen)) + baseColor;
colorMapValues = [colorMapRed colorMapGreen colorMapBlue];
colormap(colorMapValues);

if numStructures > 0
    switch hetType
        case {1, 2} % Layers or Blocks
            % Plot the X-Y plane through the center of Z.
            % Find the blocks that intersect the center of the model by
            % looking for Z values within one block length of the center.
            % If no blocks are found, move over by one gap length and try
            % again.
            structuresTemp = zeros(1,8);
            centerOfModelZ = round(numVoxelsZ/2);
            blockLengthZ = max(structures(:,7)-structures(:,4))+1;
            while sum(structuresTemp) == 0
                intersectionRangeZ = (centerOfModelZ-blockLengthZ:...
                    centerOfModelZ+blockLengthZ);
                for i = 1:length(structures(:,1))
                    if ismember(structures(i,4),intersectionRangeZ)
                        structuresTemp = [structuresTemp; structures(i,:)];
                    end
                end
                if sum(structuresTemp) == 0
                    centerOfModelZ = centerOfModelZ + (round(blockLengthZ/2));
                end
            end
            
            % Remove the first row of zeros and add one to the end indices
            % so that the plot looks like what the model should look like
            % (include the zero along the axes but also add one to reach the
            % expected value for the ending index, e.g. a width of 20
            % should align with 20 with no gap between 0 and 1)
            structuresTemp = structuresTemp(2:end,:);
            structuresTemp = [structuresTemp(:,1:4)...
                (structuresTemp(:,5:7) + 1) structuresTemp(:,8)];
            
            % Add 1 to each x value because Matlab starts with 1 not 0
            % However, in order to maintain a finite voxel width, only add one to the ending x value.
            % Example:  layer in Crystallize = 5 5, layer in Matlab = 4 5
            %           layer in Crystallize = 20 30, layer in Matlab = 19 30
            structuresTemp = [2 1 1 0 numVoxelsX numVoxelsY numVoxelsZ...
                defaultRxtAmount; structuresTemp]; % Create default rxt amount layer
            
            % Use fill to plot polygons of the layers
            for i = 1:length(structuresTemp(:,1))
                x1 = structuresTemp(i,2);
                y1 = structuresTemp(i,3);
                x2 = structuresTemp(i,5);
                y2 = structuresTemp(i,6);
                colorValue = structuresTemp(i,8);
                fill([y1 y2 y2 y1],[x1 x1 x2 x2],colorValue,'LineStyle','none')
                % coordinates are specified for all corners of polygon
                % first vector contains y coordinates, second vector contains x coordinates
                % the x and y coordinates are swapped so that the plot shows horizontal layering
                hold on
            end
            
        case 3 % Ellipsoids
            % Plot the X-Y plane through the center of Z.
            % Find the heterogeneities that intersect the center of the model by
            % looking for Z values within one het length of the center.
            structuresTemp = zeros(1,8);
            centerOfModelZ = round(numVoxelsZ/2);
            hetLengthZ = structures(2,7)+1;
            while sum(structuresTemp) == 0
                intersectionRangeZ = (centerOfModelZ-hetLengthZ:...
                    centerOfModelZ+hetLengthZ);
                for i = 1:length(structures(:,1))
                    if ismember(structures(i,4),intersectionRangeZ)
                        structuresTemp = [structuresTemp; structures(i,:)];
                    end
                end
                if sum(structuresTemp) == 0
                    centerOfModelZ = centerOfModelZ + (round(hetLengthZ/2));
                end
            end
            
            % Remove the first row of zeros and add one to the end indices
            % so that the plot looks like what the model should look like
            % (include the zero along the axes but also add one to reach the
            % expected value for the ending index, e.g. a width of 20
            % should align with 20 with no gap between 0 and 1)
            structuresTemp = structuresTemp(2:end,:);
            structuresTemp = [structuresTemp(:,1:4)...
                (structuresTemp(:,5:7) + 1) structuresTemp(:,8)];
            
            % Add 1 to each x value because Matlab starts with 1 not 0
            % However, in order to maintain a finite voxel width, only add one to the ending x value.
            % Example:  layer in Crystallize = 5 5, layer in Matlab = 4 5
            %           layer in Crystallize = 20 30, layer in Matlab = 19 30
            
            % Plot the base color (default reactant concentration)
            fill([1 numVoxelsX numVoxelsY 1],...
                [1 1 numVoxelsX numVoxelsY],defaultRxtAmount)
            % coordinates are specified for all corners of polygon
            % first vector contains y coordinates, second vector contains x coordinates
            % the x and y coordinates are swapped so that the plot shows horizontal layering
            hold on
            
            % Use fill to plot polygons of the heterogeneities
            for n = 2:length(structuresTemp(:,1))
                centerX = structuresTemp(n,2);
                radiusX = structuresTemp(n,5);
                centerY = structuresTemp(n,3);
                radiusY = structuresTemp(n,6);
                minIndexX = centerX - radiusX;
                maxIndexX = centerX + radiusX;
                minIndexY = centerY - radiusY;
                maxIndexY = centerY + radiusY;
                colorValue = structuresTemp(n,8);
                
                %check to see if any dimension is outside of the model
                %space.  If so, make the new limit the sames as the model
                %limit
                if minIndexX < 1
                    minIndexX = 1;
                end
                if minIndexY < 1
                    minIndexY = 1;
                end
                if maxIndexX > numVoxelsX
                    maxIndexX = numVoxelsX;
                end
                if maxIndexY > numVoxelsY
                    maxIndexY = numVoxelsY;
                end
                
                % Find voxels around the boundary of the heterogeneity by
                % searching for minimum and maximum coordinates along one
                % axis that falls within the heterogeneity.
                polyXY = [0 0];
                for j = minIndexY:maxIndexY
                    polyXYTemp = [0 0];
                    for i = minIndexX:maxIndexX
                        dispX = centerX - i;
                        dispY = centerY - j;
                        displacementRatio = (dispX^2)/radiusX^2 + ...
                            (dispY^2)/radiusY^2;
                        if displacementRatio <= 1
                            polyXYTemp = [polyXYTemp; j i];
                        end
                    end
                    polyXYTemp = polyXYTemp';
                    polyXYFirst = polyXYTemp(find(polyXYTemp,2));
                    polyXYLast = polyXYTemp(find(polyXYTemp,2,'last'));
                    polyXY = [polyXY; polyXYFirst'; polyXYLast'];
                end
                
                % Rearrange coordinates so that polygon has continuous
                % points around it (otherwise the outline will zigzag
                % across the polygon area).  Coordinates are acquired
                % above in pairs, so this separates the pairs into even
                % and odd sets, and then reverses the even set.
                polyXYOdd = zeros(1,2);
                polyXYEven = zeros(1,2);
                for m = 1:length(polyXY)
                    if rem(m,2) % test for odd number
                        polyXYOdd = [polyXYOdd; polyXY(m,:)];
                    else % other numbers must be even
                        polyXYEven = [polyXYEven; polyXY(m,:)];
                    end
                end
                % Reverse the order of the even array values
                polyXYEvenTemp = zeros(1,2);
                lastEven = length(polyXYEven(:,1));
                for m = lastEven:-1:2
                    polyXYEvenTemp = [polyXYEvenTemp; polyXYEven(m,:)];
                end
                polyXYEven = polyXYEvenTemp(2:end,:);
                % Put the arrays back together and remove excess zeros
                polyXY = [polyXYOdd(3:end,:); polyXYEven];
                
                %                     % Choose whether or not to outline the polygons
                %                     if outlineHet
                %                         fill(polyXY(:,1),polyXY(:,2),colorValue)
                %                     else
                fill(polyXY(:,1),polyXY(:,2),colorValue,'LineStyle','none')
                %                     end
                %                     % coordinates are specified for all corners of polygon
                %                     % first vector contains y coordinates, second vector contains x coordinates
                %                     % the x and y coordinates are swapped so that the plot shows horizontal layering
                hold on
            end
            
    end
    hold off
    set(h_rxt_axes,...
        'Position',rxtPositionSim2,...
        'XLim',[1,numVoxelsY],...
        'YLim',[1,numVoxelsX],...
        'TickDir','out',...
        'CLim',[0,4],... %Color axis limits (or Z limits)
        'FontSize',fontSizeRxtDistribution)
    title('Reactant Dist. (g/cm^3)','FontSize',...
        fontSizeRxtDistribution+2,'FontWeight','Bold');
    xlabel('y voxels')
    ylabel('x voxels')
    axis square
    colorbar
    
else % For homogeneous distributions
    fill([0 numVoxelsY numVoxelsY 0],[0 0 numVoxelsX numVoxelsX],defaultRxtAmount,'LineStyle','none')
    set(h_rxt_axes,...
        'Position',rxtPositionSim2,...
        'XLim',[1,numVoxelsY],...
        'YLim',[1,numVoxelsX],...
        'TickDir','out',...
        'CLim',[0,4],... %Color axis limits (or Z limits)
        'FontSize',fontSizeRxtDistribution)
    title('Reactant Dist. (g/cm^3)','FontSize',...
        fontSizeRxtDistribution+2,'FontWeight','Bold');
    xlabel('y voxels')
    ylabel('x voxels')
    axis square
    colorbar
    
end

% Add Crystal Size Distribution
h_CSD_axes = axes;
if rockPresent == 1  % If the natural sample is present, plot the CSD blue line for comparison with the simulation
    plot(bins,rockElements,'b','LineWidth',lineWidthForCSD);
    hold on
end
bar(bins,simulationElements,barWidthForCSD,'y')
set(h_CSD_axes,...
    'Position',csdPositionSim2,...
    'FontSize',fontSizeForCSD,...
    'TickDir','out',...
    'TickLength',lengthOfTicks,...
    'XTick',csdXTick,...
    'XLim',csdXLim,...
    'YLim',csdYLim);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
title('Crystal-Size Distribution','FontSize',fontSizeForCSD+2,'FontWeight','Bold')
if rockPresent == 1
    legend('Rock','Sim','Location','Best')
    legend('boxoff')
end
xlabel('Radius (cm)')
ylabel('Crystals (%)')
hold off

% Add nucleation rate, heating path, affinity, and Al in system
% Both max (theory) and rock-wide (simulated) nucleation rate
bothNuclPosition = ...
    [plotStackPositionSim2(1) ...
    (plotStackPositionSim2(2)+2.25*plotStackDimensionsSim2(2)) ...
    plotStackDimensionsSim2(1) plotStackDimensionsSim2(2)];
timeMyShifted = timeMy;
timeMyShifted(2:end) =  timeMy(2:end) - (reportInterval/(2*1000000));
% Shift the values by half of a reporting interval
% so that the values are plotted midway between
% reporting intervals (the number of crystals per
% second is essentially an average over the time
% interval).
h_Xtls_both_axes = axes;
h_Xtls_both = plot(timeMy,nucleationRateMax,timeMyShifted,crystalsPerSecond);
set(h_Xtls_both(1),...
    'LineStyle','-',...
    'Color','k',...
    'LineWidth',lineWidthForStack);
set(h_Xtls_both(2),...
    'LineStyle','-',...
    'Color',[0.565 0.153 0.557],... %purple
    'LineWidth',lineWidthForStack);
set(h_Xtls_both_axes,...
    'Position',bothNuclPosition,...
    'FontSize',fontSizeForStack,...
    'TickDir','out',...
    'XLim',stackXLim,...
    'XTickLabel',[],...
    'TickLength',lengthOfTicks);
ylabel('dN_V / dt')
if nucleationRateYAxis == 1
    nucleationRateYTick = nucleationRateYTick(2:end); % Remove first label to avoid
    % overlap with plot below.
    % Also add padding to make
    % space for smoothing label.
    set(h_Xtls_both_axes,...
        'YLim',nucleationRateYLim,...
        'YTick',nucleationRateYTick)
else
    nucleationRateYTick = get(h_Xtls_both_axes,'YTick');
    nucleationRateYTick = nucleationRateYTick(2:end); % Remove first label to avoid
    % overlap with plot below.
    set(h_Xtls_both_axes,'YTick',nucleationRateYTick)
end
bothTextXLim = get(h_Xtls_both_axes,'XLim');
bothTextYLim = get(h_Xtls_both_axes,'YLim');
text(bothTextXLim(2)*0.65,...
    bothTextYLim(2)*0.86,...
    'Max w/rxts',...
    'FontSize',fontSizeForStack-2)
text(bothTextXLim(2)*0.65,...
    bothTextYLim(2)*0.93,...
    'Rock-wide',...
    'FontSize',fontSizeForStack-2,...
    'Color',[0.565 0.153 0.557]) %purple
if smoothingSwitch
    text(bothTextXLim(2)*0.45,...
        bothTextYLim(2)*1.07,...
        sprintf('(%i-point smoothing)',smoothingWindowNucl),...
        'Interpreter','none','FontSize',fontSizeForStack-2,...
        'Color',[0.565 0.153 0.557]) %purple
end

% Reaction Affinity
affinityPosition = ...
    [plotStackPositionSim2(1) ...
    (plotStackPositionSim2(2)+1.25*plotStackDimensionsSim2(2)) ...
    plotStackDimensionsSim2(1) plotStackDimensionsSim2(2)];

timeMyNucl = simulationValuesXl(:,1)/10^6;

h_Affinity_axes = axes;
h_Affinity = plot(timeMyNucl,affinityNucl,...
    timeMy,affinityMean,affinityNuclFirstTime,affinityNuclFirst);
set(h_Affinity(1),...
    'LineStyle','none',...
    'Marker','o',...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor','b',...
    'MarkerSize',2);
set(h_Affinity(2),...
    'LineStyle','-',...
    'Color','k',...
    'LineWidth',1);
set(h_Affinity(3),...
    'LineStyle','none',...
    'Marker','o',...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','none',...
    'MarkerSize',2);
set(h_Affinity_axes,...
    'Position',affinityPosition,...
    'FontSize',fontSizeForStack,...
    'TickDir','out',...
    'XLim',stackXLim,...
    'XTickLabel',[],...
    'TickLength',lengthOfTicks);
ylabel('-\Delta_rG')
if affinityYAxis == 1
    affinityYTick = affinityYTick(2:end);
    set(h_Affinity_axes,...
        'YLim',affinityYLim,...
        'YTick',affinityYTick)
else
    affinityYTick = get(h_Affinity_axes,'YTick');
    affinityYTick = affinityYTick(2:end);
    set(h_Affinity_axes,...
        'YLim',affinityYLim,...
        'YTick',affinityYTick)
end

% Heating Path
heatingPathPosition = ...
    [plotStackPositionSim2(1) ...
    (plotStackPositionSim2(2)+0.75*plotStackDimensionsSim2(2)) ...
    plotStackDimensionsSim2(1) plotStackDimensionsSim2(2)*0.50];
h_HeatPath_axes = axes;
h_HeatPath = plot(timeMy,temperature);
set(h_HeatPath,...
    'LineStyle','-',...
    'Color','k',...
    'LineWidth',1)
set(h_HeatPath_axes,...
    'Position',heatingPathPosition,...
    'FontSize',fontSizeForStack,...
    'TickDir','out',...
    'XLim',stackXLim,...
    'YLim',heatingPathYlim,...
    'YTick',heatingPathYTick,...
    'XTickLabel',[],...
    'TickLength',lengthOfTicks);
ylabel('T (\circC)')

% Percent Al in CAR and CAP
systemAlPosition = ...
    [plotStackPositionSim2(1) ...
    (plotStackPositionSim2(2)+0.25*plotStackDimensionsSim2(2)) ...
    plotStackDimensionsSim2(1) plotStackDimensionsSim2(2)*0.50];
h_Al_axes = axes;
h_Al = plot(timeMy,reactantAl,timeMy,productAl);
set(h_Al(1),...
    'LineStyle','-',...
    'Color','g',...
    'LineWidth',lineWidthForStack)
set(h_Al(2),...
    'LineStyle','-',...
    'Color','r',...
    'LineWidth',lineWidthForStack)
set(h_Al_axes,...
    'Position',systemAlPosition,...
    'FontSize',fontSizeForStack,...
    'TickDir','out',...
    'XLim',stackXLim,...
    'XTickLabel',[],...
    'YLim',percentAlYLim,...
    'YTick',percentAlYTick,...
    'TickLength',lengthOfTicks);
ylabel('Al (%)')

% Percent Al in Fluid
h_Al_fluid_axes = axes;
h_Al_fluid = plot(timeMy,fluidAl);
set(h_Al_fluid,...
    'LineStyle','-',...
    'Color','cyan',...
    'LineWidth',lineWidthForStack)
set(h_Al_fluid_axes,...
    'Position',[plotStackPositionSim2 plotStackDimensionsSim2(1)...
    plotStackDimensionsSim2(2)*0.25],...
    'FontSize',fontSizeForStack,...
    'TickDir','out',...
    'XLim',stackXLim,...
    'TickLength',lengthOfTicks);
xlabel('Time (Myr)')
ylabel('Al (%)')



% Write the simulation title, simulation parameters, and rock parameters

% Create axes for entire figure (need this to write text to fig)
axes('position',[0 0 1 1],'visible','off')

% Title
if rockPresent == 1
    text(titleRockPositionSim2(1),titleRockPositionSim2(2),rockName,...
        'FontWeight','Bold','FontSize',16)
end
text(titleSimPositionSim2(1),titleSimPositionSim2(2),...
    simulationName,'Interpreter','none',...
    'FontWeight','Bold','FontSize',12)

% Calculate the text positions
yPosition = textPositionSim2(2):...
    textIncrementSim2:...
    textPositionSim2(2)+(numberOfLinesSim2-1)*textIncrementSim2;


% Descriptions
text(xDescPositionSim2,yPosition(28)+textGapSim2*3,'N, total crystals:'                      ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(27)+textGapSim2*3,'N_{max}, specified maximum number of xtls:','FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(26)+textGapSim2*3,'ND, xtl num density (xtls/cm^3):'        ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(25)+textGapSim2*3,'t_{dur}, xtlzn duration (my):'           ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(24)+textGapSim2*3,'T_{eq}, equil T of rxn (\circC):'        ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(23)+textGapSim2*3,'T_{over}, therm overstep (\circC):'      ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(22)+textGapSim2*3,'T_{95}, T at 95% Al in prod (\circC):'   ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(21)+textGapSim2*3,'r_{min}, min radius (cm):'               ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(20)+textGapSim2*3,'r_{mean}, mean radius (cm):'             ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(19)+textGapSim2*3,'r_{max}, max radius (cm):'               ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(18)+textGapSim2*3,'Mode (vol%):'                            ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(17)+textGapSim2*3,'Max mode, specified (ppb vol%):'         ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(16)+textGapSim2*2,'D_{inf} (cm^2/s):'                       ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(15)+textGapSim2*2,'[Al]_{fl} in eq w/products (mol/cm^3):'  ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(14)+textGapSim2*2,'\phi, porosity:'                         ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(13)+textGapSim2*2,'\tau, tortuosity:'                       ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(12)+textGapSim2*2,'Q_D (kJ/mol):'                           ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(11)+textGapSim2*2,'D = D_{inf} \phi \tau e^{(-Q_D/RT)} (m^2/s):','FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(9)+textGapSim2*2,'k_1, (dN/dt)_{steady-state} (nucl/cm^3/s):','FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(8)+textGapSim2*2,'(dN/dt)_{max} (nucl/cm^3/s):'             ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(7)+textGapSim2*2,'k_2, nucl acceleration:'                  ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(6)+textGapSim2,'A_1, first nucl affinity (kJ/mol):'         ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(5)+textGapSim2,'A_{max}, max nucl affinity (kJ/mol):'       ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(4)+textGapSim2,'A_{max mean}, max mean affinity (kJ/mol):'  ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(3),'Time step (y):'                                     ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(2),'Reporting interval (y):'                            ,'FontSize',sizeOfText)
text(xDescPositionSim2,yPosition(1),'^aDetermined at nearest reporting interval'         ,'FontSize',sizeOfText-1)

% Model and simulation parameters
text(xSimPositionSim2,yPosition(28)+textGapSim2*3,sprintf('%0.0f',length(simulationValuesInt)),'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(27)+textGapSim2*3,sprintf('%0.0f',max_ppb)                   ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(26)+textGapSim2*3,sprintf('%0.0f',crystalNumberDensitySim)   ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(25)+textGapSim2*3,sprintf('%2.1f^a',timeRange)               ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(24)+textGapSim2*3,sprintf('%3.0f',Teq)                       ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(23)+textGapSim2*3,sprintf('%3.0f',delaT)                     ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(22)+textGapSim2*3,sprintf('%3.0f^a',temperature95)           ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(21)+textGapSim2*3,sprintf('%5.3f',simulationRadiiMin)        ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(20)+textGapSim2*3,sprintf('%5.3f',simulationRadiiMean)       ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(19)+textGapSim2*3,sprintf('%5.3f',simulationRadiiMax)        ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(18)+textGapSim2*3,sprintf('%3.1f',simulationPorphMode)       ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(17)+textGapSim2*3,sprintf('%3.1f',max_mode_ppb)              ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(16)+textGapSim2*2,sprintf('%4.2e',Dinf)                      ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(15)+textGapSim2*2,sprintf('%4.2e',CiAl)                      ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(14)+textGapSim2*2,sprintf('%4.2e',porosity)                  ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(13)+textGapSim2*2,sprintf('%4.2e',tortuosity)                ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(12)+textGapSim2*2,sprintf('%g',Q)                            ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(11)+textGapSim2*2,sprintf('%4.2e^a',D)                       ,'FontName','FixedWidth')
text(xDescPositionSim2+0.020,yPosition(10)+textGapSim2*2,sprintf('(at T_c, %0.0f ^oC and %3.1f Myr)^a',tempD,timeD),'FontSize',sizeOfText)
text(xSimPositionSim2,yPosition(9)+textGapSim2*2,sprintf('%4.2e',k1)                         ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(8)+textGapSim2*2,sprintf('%4.2e',max(nucleationRateMax))     ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(7)+textGapSim2*2,sprintf('%4.2e',k2)                         ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(6)+textGapSim2,sprintf('%2.1f',affinityNuclFirst)            ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(5)+textGapSim2,sprintf('%2.1f',affinityNuclMax)              ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(4)+textGapSim2,sprintf('%2.1f',affinityMeanMax)              ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(3),sprintf('%g',timeStep)                                ,'FontName','FixedWidth')
text(xSimPositionSim2,yPosition(2),sprintf('%4.2e',reportInterval)                       ,'FontName','FixedWidth')

% Rock parameters
if rockPresent == 1
    text(xRockPositionSim2,yPosition(numberOfLinesSim2)+textGapSim2*5                      ,'Rock','FontWeight','Bold')
    text(xSimPositionSim2,yPosition(numberOfLinesSim2)+textGapSim2*5                 ,'Simulation','FontWeight','Bold')
    text(xRockPositionSim2,yPosition(28)+textGapSim2*3,num2str(length(rockValuesInt(:,1))),'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(26)+textGapSim2*3,crystalNumberDensityRock           ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(25)+textGapSim2*3,timeRangeString                    ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(24)+textGapSim2*3,temperatureEquilibriumRock         ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(22)+textGapSim2*3,temperature95Rock                  ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(21)+textGapSim2*3,sprintf('%5.3f',rockRadiiMin)      ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(20)+textGapSim2*3,sprintf('%5.3f',rockRadiiMean)     ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(19)+textGapSim2*3,sprintf('%5.3f',rockRadiiMax)      ,'FontName','FixedWidth')
    text(xRockPositionSim2,yPosition(18)+textGapSim2*3,rockMode                           ,'FontName','FixedWidth')
end


% Make a 3D plot of the crystals
if plotCrystalsIn3D
    figure('FileName',strcat(simulationName,'_Plot3'))
    % Use bubbleplot3 function to create a plot
    bubbleValuesInt = simulationValuesInt';
    bubbleplot3(bubbleValuesInt(2,:),bubbleValuesInt(3,:),bubbleValuesInt(4,:),bubbleValuesInt(5,:));
    % Set the figure properties
    title(simulationName,'Interpreter','none','FontWeight','Bold','FontSize',12)
    xMax = ceil((voxelEdgeLength*numVoxelsX));
    yMax = ceil((voxelEdgeLength*numVoxelsY));
    xlim([0 (voxelEdgeLength*numVoxelsX)])
    ylim([0 (voxelEdgeLength*numVoxelsY)])
    set(gca,'xtick',(0:xMax/10:xMax),'ytick',(0:yMax/10:yMax));
    xlabel('X (cm)')
    ylabel('Y (cm)')
    zlabel('Z (cm)')
    view(90,-90) % (azimuth,elevation from X-Y plane to Z plane.  90 = looking along Z)
    axis square
    %axis tight
    grid on
    camlight left
    lighting gouraud %phong
    %shading interp
    % Add a colorbar to show depth
    h_colorbar_axes = colorbar;
    set(get(h_colorbar_axes,'YLabel'),'String','Z (cm)')
    colorMapValues = zeros(181,3);
    colorMapValues(:,1) = (1:-0.005:0.1);
    colorMapValues(:,2) = 0.1;
    colorMapValues(:,3) = 0.1;
    colormap(colorMapValues);
    
end
% Plot zoning vs radius relationship (Sim4 is the identifier)
if plotZoning
    % Figure position
    figure('Position',[(figPositionSim2(1)+figPositionSim2(3)+50)...
        figPositionSim2(2) figPositionSim2(3) figPositionSim2(4)],...
        'PaperPosition',[0.25,0.25,8,10.5],...
        'FileName',strcat(simulationName,'_Plot4'))
    
    simColor = [1 0.5 0]; % orange
    
    % Core concentration vs final radius
    subplot(3,4,1:2)
    plot(rockValuesRadiusVsCoreConc(:,1),rockValuesRadiusVsCoreConc(:,2),...
        'bd','MarkerFaceColor','b','MarkerSize',7)
    hold on
    plot(simulationRadiusVsCoreConc(:,1),simulationRadiusVsCoreConc(:,2),...
        'ro','MarkerFaceColor',simColor,'MarkerSize',5,'MarkerEdgeColor','k')
    hold off
    set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
    set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
    set(gca,'FontSize',8)
    title({simulationName},'Interpreter','none','FontSize',10,...
        'FontWeight','Bold')
    xlabel('Final Radius (cm)')
    ylabel('Core Concentration (Wt%)')
    legend('Rock','Sim','Location','NorthWest')
    legend('boxoff')
    
    % Radius range vs core concentration (radius range is a proxy for
    % competition)
    subplot(3,4,3)
    semilogx(ratioOfRadiusRange,culledBins)
    title('Competition')
    grid on
    xlabel('\Deltar_{sim}/\Deltar_{rock}')
    set(gca,'XLim',[0.1 10],...
        'XTick',[0.1 1 10],...
        'XMinorGrid','off',...
        'YTick',binsCoreConc,...
        'YGrid','off',...
        'YTickLabel',[],...
        'FontSize',8);
        %'XTickLabel',num2str(get(gca,'XTick')','%.0f'));
        
    % Volume percent vs core concentration (volume is a proxy for growth
    % rate and additional transport mechanisms, like fluid flux)
    subplot(3,4,4)
    semilogx(ratioOfVol,culledBins);
    title('Growth')
    grid on
    xlabel('V%_{sim/rock}')
    set(gca,'XLim',[0.1 10],...
        'XTick',[0.1 1 10],...
        'XMinorGrid','off',...
        'YTick',binsCoreConc,...
        'YGrid','off',...
        'YTickLabel',[],...
        'FontSize',8);
         %'XTickLabel',num2str(get(gca,'XTick')','%.0f'),...
   
    % Histogram of core concentration
    subplot(3,4,5:6)
    elementsRock = hist(rockValuesRadiusVsCoreConc(:,2),binsCoreConc);
    elementsRock = 100*(elementsRock./sum(elementsRock));
    elementsSim = hist(simulationRadiusVsCoreConc(:,2),binsCoreConc);
    elementsSim = 100*(elementsSim./sum(elementsSim));
    plot(binsCoreConc,elementsRock,'b-','LineWidth',2)
    hold on
    bar(binsCoreConc,elementsSim,0.7,'FaceColor',simColor)
    hold off
    xlim(zoningVsRadiusXLim)
    set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.1f'),'FontSize',8)
    title('Relative nucleation time','FontWeight','Bold')
    xlabel('Core Concentration (Wt%)')
    ylabel('Crystals (%)')
    legend('Rock','Sim','Location','NorthEast')
    legend('boxoff')
       
    % Plot the volumetric growth rate of the first, middle, and last
    % crystals
    % Determine mid-points between time steps and calculate differences
    % between time steps (last step is sometimes not the same)
    simulationTimeMidPoints = zeros(1,length(simulationValuesRadiiVsTime(1,:)));
    simulationTimeDiffs = zeros(size(simulationTimeMidPoints));
    for i = 2:(length(simulationValuesRadiiVsTime(1,:)))
        simulationTimeMidPoints(i) = (simulationValuesRadiiVsTime(1,i) + ...
            simulationValuesRadiiVsTime(1,i-1)) / 2;
        simulationTimeDiffs(i) = (simulationValuesRadiiVsTime(1,i) - ...
            simulationValuesRadiiVsTime(1,i-1));
    end
    % Remove first column of excess zeros
    simulationTimeMidPoints = simulationTimeMidPoints(1,2:end);
    simulationTimeDiffs = simulationTimeDiffs(1,2:end);
    % Convert radii to volumes and calculate the volumetric growth rates
    simulationValuesVol = ((4/3)*pi).*(simulationValuesRadiiVsTime(2:end,:).^3);
    simulationValuesVolVsTime = [simulationValuesRadiiVsTime(1,:); simulationValuesVol];
    simulationValuesVRateVsTime = zeros(size(simulationValuesVolVsTime));
    for i = 2:length(simulationValuesVolVsTime(:,1))
        for j = 2:length(simulationValuesVolVsTime(1,:))
            simulationValuesVRateVsTime(i,j) = ...
                (simulationValuesVolVsTime(i,j) - simulationValuesVolVsTime(i,j-1)) /...
                simulationTimeDiffs(j-1);
        end
    end
    % Remove first row and column of excess zeros
    simulationValuesVRateVsTime = simulationValuesVRateVsTime(2:end,2:end);
    simulationValuesVRateVsTime = [simulationTimeMidPoints; ...
        simulationValuesVRateVsTime];
    % Determine the index of the middle crystal
    middleCrystalIndex = round(length(simulationValuesVRateVsTime(:,1)) / 2);
    % Plot the values
    subplot(3,4,7:8)
    plot(simulationValuesVRateVsTime(1,:),simulationValuesVRateVsTime(2,:),...
        simulationValuesVRateVsTime(1,:),simulationValuesVRateVsTime(middleCrystalIndex,:),...
        simulationValuesVRateVsTime(1,:),simulationValuesVRateVsTime(end,:))
    title('Simulated growth rate','FontSize',8,'FontWeight','Bold')
    legend('First','Middle','Last','Location','NorthWest')
    legend('boxoff')
    set(gca,'FontSize',8)
    xlabel('Time (Myr)','FontSize',8)
    ylabel('Volumetric rate (cm^3/Myr)','FontSize',8)

    % For the spline plots, calculate the values of the time vs zoning fit
    % to the radii-vs-time simulation output.
    simulationZoningFromTime = ...
        ppval(simulationTimeAndZoningSplineCoef,simulationTimeAndZoning(:,1));
    simulationZoningFromTimeSmoothed = ...
        smooth(simulationZoningFromTime,smoothingWindowZoning);
    
    % Zoning as a function of radius for the natural Grt and the spline
    % expressed in two plots
    subplot(3,4,9:10)
    plot(rockValuesGrtZoningNorm(:,1),rockValuesGrtZoningNorm(:,2),'bo',...
        'MarkerSize',3)
    hold on
    plot(simulationFirstRadiiVsTime(:,2),rockZoningFromRadiusSmoothed,'r-')
    hold off
    garnetZoningTitle = 'Measured garnet zoning profile';
    title(garnetZoningTitle,...
        'Interpreter','none','FontSize',8,'FontWeight','Bold')
    set(gca,'FontSize',8)
    xlabel('Radius (cm)','FontSize',8)
    ylabel('Concentration','FontSize',8)
    if smoothingWindowZoning > 0
        h_spline1 = legend('Measured',...
            sprintf('Spline fit (%i-pt smoothing)',smoothingWindowZoning));
    else
        h_spline1 = legend('Measured','Spline fit');
    end
    set(h_spline1,'Location','SouthWest','Box','off','FontSize',7)
    
    % Rock-wide concentration
    subplot(3,4,11:12)
    plot(simulationTimeAndZoning(:,1),simulationTimeAndZoning(:,2),'bo',...
        'MarkerSize',3)
    hold on
    plot(simulationTimeAndZoning(:,1),simulationZoningFromTimeSmoothed,'r-')
    hold off
    title('Calculated rock-wide concentration',...
        'FontSize',8,'FontWeight','Bold')
    set(gca,'FontSize',8)
    xlabel('Time (Myr)','FontSize',8)
    ylabel('Concentration','FontSize',8)
    if smoothingWindowZoning > 0
        h_spline2 = legend('Simulation',...
            sprintf('Spline fit (%i-pt smoothing)',smoothingWindowZoning));
    else
        h_spline2 = legend('Simulation','Spline fit');
    end
    set(h_spline2,'Location','SouthWest','Box','off','FontSize',7)
    
    % Attach the underlying data to the models in the subplot
    axes('Position',[0 0 1 1],'Visible','off');
    text(0.025,0.050,sprintf('Crystal sizes: %s',rockFileName),...
        'Interpreter','none','FontSize',7)
    text(0.025,0.035,sprintf('Core-rim Grt zoning file: %s',grtZoningFileName),...
        'Interpreter','none','FontSize',7)
    text(0.025,0.020,sprintf('Core concentrations and radii: %s',radiusVsCoreConcFileName),...
        'Interpreter','none','FontSize',7)
end