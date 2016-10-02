function [layerBoundariesVoxels,crystalRadius] =...
    calcCrystalsInLayers(numVoxels,sizeVoxel,simValues,structures)

% calcCrystalsInLayers.m
% Function for extracting the crystals in each layer
% The crystals in each layer can be plotted as CSDs and
% used to understand the contribution of each layer to
% the whole-rock CSD.  This is useful in the iterative
% process of determining the initial reactant
% distribution for a model.


% Extract layer boundaries as voxel numbers
%layerBoundariesVoxels = structures(:,2:3);

% Find and add the non-explicit layers

% sort layers by column 1
layerBoundariesVoxels = sortrows(structures(:,2:3),1);

if layerBoundariesVoxels(1,1) > 0
    % make temp matrix starting at 0
    layerBoundariesVoxelsTemp = [0 (layerBoundariesVoxels(1,1)-1)];
    % append current layer to the temp matrix
    layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;layerBoundariesVoxels(1,:)];
else
    % make temp matrix using first row
    layerBoundariesVoxelsTemp = layerBoundariesVoxels(1,:);
end

if length(structures(:,1)) > 1
	for i=2:length(structures(:,1))
        if (layerBoundariesVoxels(i,1) - layerBoundariesVoxels(i-1,2)) == 1
            % append current layer to the temp matrix
            layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;layerBoundariesVoxels(i,:)];
        elseif layerBoundariesVoxels(i,1) < layerBoundariesVoxels(i-1,2)
            if layerBoundariesVoxels(i,1) == 0
                % append current layer to the temp matrix
                layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;layerBoundariesVoxels(i,:)];
            else
                % append new layer to the temp matrix
                layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;...
                        0 layerBoundariesVoxels(i,1)-1];
                % append current layer to the temp matrix
                layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;layerBoundariesVoxels(i,:)];
            end
        else
            % append new layer to the temp matrix
            layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;...
                    layerBoundariesVoxels(i-1,2)+1 layerBoundariesVoxels(i,1)-1];
            % append current layer to the temp matrix
            layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;layerBoundariesVoxels(i,:)];
        end
	end
end

if layerBoundariesVoxels(end,2) < (numVoxels-1)
    % append new layer to the temp matrix
    layerBoundariesVoxelsTemp = [layerBoundariesVoxelsTemp;...
            layerBoundariesVoxels(end,2)+1 numVoxels-1];
end

layerBoundariesVoxels = layerBoundariesVoxelsTemp;

% Find the crystals that fall within each layer
% This compares the x coordinate of the crystal with the 
% x-coordinate boundaries of the layers
layerBoundariesCM = sizeVoxel.*layerBoundariesVoxels; % Convert to centimeters
crystalRadius = zeros(length(simValues(:,1)),1); % Establish a matrix to hold the values
radius = crystalRadius; % Establish a matrix to hold the values
for i=1:length(layerBoundariesCM(:,1))
    for j=1:length(simValues(:,1))
        if (layerBoundariesCM(i,1) <= simValues(j,2)) && (simValues(j,2) <= layerBoundariesCM(i,2))
           radius(j) = simValues(j,5);
       else
           radius(j) = 0;
        end
    end
    crystalRadius = [crystalRadius radius]; % Append each new column of radii to the matrix
end
crystalRadius = crystalRadius(:,2:end); % Remove the first column of zeros created above

