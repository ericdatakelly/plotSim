function [density] = calcCrystalNumDensity(values,volume);

% calcCrystalNumDensity.m
% This funciton calculates the number of crystals per volume

crystalsTotal = max(values(:,10));
density = crystalsTotal/volume;