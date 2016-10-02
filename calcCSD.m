function [elements] = calcCSD(values,bins)

% calcCSD.m
% This function extracts and normallizes the values used in the crystal-size distribution plot
% Elements are the values in a histogram and the bins are the categories that the elements fit into

elements = hist(values,bins);
elements = 100*(elements./sum(elements));
