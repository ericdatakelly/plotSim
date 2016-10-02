function [timeMy,crystals] = ...
    calcCrystalsPerTime(volume,values,smoothingSwitch,windowSize)

% calcCrystalsPerTime.m
% This function sums the crystals at each reporting interval
% and divides by the volume of the sample to produce
% the number of crystals per cm^3 per interval
% Time is converted into millions of years


% Calculate Number of Crystals at each time interval
	timeMy = values(:,2)/1000000;
    crystals = zeros(size(values(:,10)));
    crystals(1) = values(1,10);
    for i = 2:length(values(:,10))
        crystals(i) = (values(i,10) - values(i-1,10));
    end
    timeStep = timeMy(2) - timeMy(1);
    crystals = crystals(2:end)./(timeStep * 3.1557e13 * volume); %Normalize to seconds and cm^3
    crystals = [0;crystals]; % Pad with 0 to maintain length and avoid 
                             % mismatched arrays in other functions.
    
    % The following will smooth the values based on the user-defined window
    % in the plotting parameters file
    if smoothingSwitch
        %crystals = filter(ones(1,windowSize)/windowSize,1,crystals)
        crystals = smooth(crystals,windowSize);
    end

end