function [binnedRadiusRangeAndVolPerc] = ...
    calcRadRangeVol(radiusConcValues,binsOfConc,vol)

% calcRadRangeVol.m
% This function extracts the minimum and maximum radius within specified
% ranges of core concentration (bins) and also determines the volume percent
% of the crystals within the binned range.
    
% Determine the volume of the Rock and Sim crystals in each bin
binnedRadiusRangeAndVol = [0 0 0 0];
% Step through all bins to find crystals that fit
for i = 1:length(binsOfConc)-1
    minRadius = max(radiusConcValues(:,1));
    maxRadius = 0;
    volTemp = 0;
    % Scan all crystals within the bin
    for j = 1:length(radiusConcValues(:,2))
        % Determine whether the crystal concentration fits the bin
        if (binsOfConc(i) <= radiusConcValues(j,2)) && ...
                (radiusConcValues(j,2) < binsOfConc(i+1))
            % Look for the smallest crystal in the bin
            if radiusConcValues(j,1) < minRadius
                minRadius = radiusConcValues(j,1);
            end
            % Look for the largest crystal in the bin
            if radiusConcValues(j,1) > maxRadius
                maxRadius = radiusConcValues(j,1);
            end
            % Sum the volumes of all crystals in the bin
            volTemp = volTemp + ((4/3) * pi * radiusConcValues(j,1)^3);
        end
    end
    % Record the values from each bin in rows of an array.  Use the average
    % value between the bin values so that the radii and volumes are
    % associated with the center of the bin.
    binnedRadiusRangeAndVol = [binnedRadiusRangeAndVol; ...
        ((binsOfConc(i)+binsOfConc(i+1))/2) minRadius maxRadius volTemp];
end

% Remove the first row of zeros and remove bins with no crystals
binnedRadiusRangeAndVol = binnedRadiusRangeAndVol(2:end,:);
% goodBins = [];
% for i = 1:length(binsOfConc)
%     if binnedRadiusRangeAndVol(i,4) > 0
%         goodBins = [goodBins i];
%     end
% end
% binnedRadiusRangeAndVol = binnedRadiusRangeAndVol(goodBins,:);

% Determine volume percent of crystals in rock or simulation
binnedRadiusRangeAndVolPerc = [binnedRadiusRangeAndVol(:,1:3) 100 *...
    (binnedRadiusRangeAndVol(:,4) / vol)];
end