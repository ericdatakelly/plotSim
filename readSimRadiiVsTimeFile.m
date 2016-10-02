function [values] = readSimRadiiVsTimeFile(name)

% readSimRadiiVsTimeFile.m
% This function reads the values from a XllizeViz radii-vs-time file
% 'name' is the path and file name of the file

fid = fopen(name);
values = readtable(name,'Delimiter','\t','HeaderLines',2,...
    'ReadVariableNames',false);
values = values{:,5:end};
values = (values')'; % This helps to convert the table to an array that I can use.
fclose(fid);