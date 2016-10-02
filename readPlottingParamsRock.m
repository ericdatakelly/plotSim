function [rockFileName,...
        grtZoningFileName,...
        radiusVsCoreConcFileName,...
        crystalNumberDensityRock,...
        timeRangeString,...
        temperatureEquilibriumRock,...
        temperature95Rock,...
        rockMode,...
        rockVol,...
        binMin,...
        binCenter,...
        binMax,...
        csdXLimMin,...
        csdXLimMax,...
        csdXTickMin,...
        csdXTickIncrement,...
        csdXTickMax,...
        csdYLimMin,...
        csdYLimMax,...
        plotLayerCSDs,...
        stackXLimMin,...
        stackXLimMax,...
        nucleationRateYAxis,...
        nucleationRateYLimMin,...
        nucleationRateYTickInc,...
        nucleationRateYLimMax,...
        smoothingSwitch,...
        smoothingWindowNucl,...
        heatingPathYlimMin,...
        heatingPathYlimInc,...
        heatingPathYlimMax,...
        plotCrystalsIn3D,...
        affinityYAxis,...
        affinityYLimMin,...
        affinityYLimInc,...
        affinityYLimMax,...
        plotZoning,...
        coreConcXMin,...
        coreConcXMax,...
        coreConcTickXMin,...
        coreConcTickXInc,...
        coreConcTickXMax,...
        smoothingWindowZoning...
] = readPlottingParamsRock(name)

% readPlottingParamsRock.m
% This function reads the plotting parameters file for the given rock
% These variables are different types, e.g. strings and numbers, so I did not store the
% values in a matrix.  There's a way to do this but I won't bother right now.
% 'name' is the path and file name of the plotting params file


fid = fopen(name);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line)); % find the colon and get everything after
rockFileName = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
grtZoningFileName = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
radiusVsCoreConcFileName = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
crystalNumberDensityRock = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
timeRangeString = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
temperatureEquilibriumRock = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
temperature95Rock = sscanf(line(2:length(line)),'%c');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
rockMode = sscanf(line(2:length(line)),'%c');line = fgetl(fid);
line = line(strfind(line,':'):length(line));
rockVol = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
binMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
binCenter = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
binMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdXLimMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdXLimMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdXTickMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdXTickIncrement = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdXTickMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdYLimMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
csdYLimMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
plotLayerCSDs = sscanf(line(2:length(line)),'%i');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
stackXLimMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
stackXLimMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
nucleationRateYAxis = sscanf(line(2:length(line)),'%i');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
nucleationRateYLimMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
nucleationRateYTickInc = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
nucleationRateYLimMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
smoothingSwitch = sscanf(line(2:length(line)),'%i');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
smoothingWindowNucl = sscanf(line(2:length(line)),'%i');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heatingPathYlimMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heatingPathYlimInc = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heatingPathYlimMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
plotCrystalsIn3D = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
affinityYAxis = sscanf(line(2:length(line)),'%i');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
affinityYLimMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
affinityYLimInc = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
affinityYLimMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
plotZoning = sscanf(line(2:length(line)),'%i');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
coreConcXMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
coreConcXMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
coreConcTickXMin = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
coreConcTickXInc = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
coreConcTickXMax = sscanf(line(2:length(line)),'%g');
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
smoothingWindowZoning = sscanf(line(2:length(line)),'%i');
fclose(fid);
