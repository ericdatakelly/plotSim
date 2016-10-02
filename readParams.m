function [reportInterval,...
    porosity,...
    tortuosity,...
    Q,...
    Dinf,...
    Teq,...
    CiAl,...
    k2,...
    k1,...
    max_ppb,...
    max_mode,...
    voxelEdgeLength,...
    timeStep,...
    numVoxelsX,...
    numVoxelsY,...
    numVoxelsZ,...
    defaultRxtAmount,...
    numStructures,...
    hetType,...
    structures,...
    fv,...
    n_GtoC...
    ] = readParams(name)

% readParams.m
% This function reads the values from the CRYSTALLIZE parameters file
% The values can probably be stored in a matrix but I haven't done it yet
% 'name' is the path and file name of the parameters file


fid = fopen(name);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line)); % find the ":" and get everything after
reportInterval = sscanf(line,'%*c %g',1);% skip the characters and get the number
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
porosity = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
tortuosity = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
Q = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
Dinf = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
CiAl = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
k2 = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
k1 = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
max_ppb = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
max_mode = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
voxelEdgeLength = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
timeStep = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
numVoxelsX = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
numVoxelsY = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
numVoxelsZ = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
defaultRxtAmount = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
numStructures = sscanf(line,'%*c %i',1);
hetType = 0;
if numStructures > 0
    structures = zeros(length(numStructures),8);
    for i = 1:numStructures
        line = fgetl(fid);
        line = line(strfind(line,')'):length(line));
        heterogeneity = sscanf(line,'%*c %g');
        switch heterogeneity(1,1)
            case 1 % 1 = layer
                hetType = 1;
                structures(i,1) = heterogeneity(1,1);
                structures(i,2) = heterogeneity(2,1);
                structures(i,5) = heterogeneity(3,1);
                structures(i,6) = numVoxelsY-1;
                structures(i,7) = numVoxelsZ-1;
                structures(i,8) = heterogeneity(4,1);
            case 2 % 2 = block
                hetType = 2;
                structures(i,:) = heterogeneity;
            case 3 % 3 = ellipse
                hetType = 3;
                structures(i,:) = heterogeneity;
        end
    end
else
    structures = zeros(1,8);
end
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
fv = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
n_GtoC = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heatingPath = sscanf(line,'%*c (%g,%g,%g)');
heatingPath = heatingPath';
Teq = heatingPath(1,2);
fclose(fid);
end