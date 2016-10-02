function [values] = readRockZoningFile(name)

% readRockZoningFile.m
% This function reads the values from two columns of data
% All values are read except for the header
% 'name' is the path and file name

fid = fopen(name);
header = str2mat(fgetl(fid));
values = fscanf(fid,'%g %g',[2 inf]);
values = values';
fclose(fid);