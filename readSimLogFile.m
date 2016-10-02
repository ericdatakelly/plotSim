function [courant] = readSimLogFile(name)

% readSimLogFile.m
% This function reads the log file from CRYSTALLIZE
% At this point it only reads the one value
% If needed, the rest of the file could be read but
% it does not seem necessary right now
% 'name' is the path and file name of the log file

fid = fopen(name);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(findstr(line,'='):length(line));
courant = sscanf(line,'%*c %g',1);

% I'm not sure how to extract the rest of the log file.  Perhaps like this...
% while findstr(line,'Time')>0
%     simulationValuesLog = sscanf(line,'Time = %g, delC = %g, NuclRate = %g',[3 inf]);
%     line = fgetl(fid);
% end
% simulationValuesLog = simulationValuesLog';
fclose(fid);