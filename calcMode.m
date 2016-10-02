function [mode] = calcMode(volume,values)

% calcMode.m
% This function calculates the mode of the crystals in the simulation
% This could be adapted to calculate the mode of the crystals in the
% rock also, but I am using the values output from REDUCE to report the
% mode of the rock so this function is not sep up to use that information


mode = 100.*sum(values)/volume;