function [values] = readSimRunFile(name)

% readSimValuesRun.m
% This function reads the values from a CRYSTALLIZE run file
% All values are read except for the header
% 'name' is the path and file name of the run file


fid = fopen(name);
header = char(fgetl(fid));
values = fscanf(fid,'%g %g %g %g %g %g %g %g %g %i %g',[11 inf]);
values = values';
fclose(fid);
end