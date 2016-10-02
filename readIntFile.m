function [values] = readIntFile(name)

% readIntFile.m
% This function reads all of the values in an integrate file
% The header information is read but not returned to the calling
% function because there is no use for it (at least for now)
% Integrate files contain the coordinates and radius values for
% a rock or a CRYSTALLIZE simulation
% 'name' is a string containing the path and file name of the int file
% The values in the file are returned as a matrix with 8 columns
% and an unlimited number of rows


fid = fopen(name);
header = fgetl(fid);
for n=1:4
    header = str2mat(header,fgetl(fid));
end
values = fscanf(fid,'%i %g %g %g %g %g %g',[7 inf]);
values = values';
fclose(fid);
