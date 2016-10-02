function [totalVolume] = calcVolume(edgeLength,x,y,z)

% calcVolume.m
% This function calculates the volume of the simulation space


voxelVolume = edgeLength^3;
totalVoxels = x*y*z;
totalVolume = voxelVolume*totalVoxels;
