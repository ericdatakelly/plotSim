function [values] = readSimXlFile(name)

% readSimValuesRun.m
% This function reads the values from a CRYSTALLIZE xl file
% All values are read except for the header
% 'name' is the path and file name of the run file

fid = fopen(name);

header = str2mat(fgetl(fid));
header = str2mat(fgetl(fid));
values = fscanf(fid,'%g %g %g %g %g %g %g %g',[8 inf]);
values = values';

fclose(fid);