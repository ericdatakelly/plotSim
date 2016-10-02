function [simulationRadiusAndCoreConc,rockValuesGrtZoningNorm,...
    rockZoningFromRadiusSmoothed,simulationTimeAndZoning,...
    simulationTimeAndZoningSplineCoef] = ...
    calcZoningVsRadius(simulationValuesXlTimeVsSize,rockValuesGrtZoning,...
    simulationValuesTimeVsRadius,smoothingWindowZoning)

% calcMnVsRadius.m
% Function for plotting chemical zoning versus radius relationships. 
% Currently, MnO is the zoning of interest.
% The basic idea is to use a measure core-rim zoning profile from a large garnet
% with high core MnO, normalize the profile to the radius of the first 
% crystal in the simulation, associate time with each radius of the
% simulated crystal (this also associates zoning with time in the simulation),
% and then use the nucleation time of each crystal in the model to assign a
% core zoning concentration.  In plotSim.m, the core zoning concentrations with their
% corresponding final radii are plotted on top of the measured core-zoning
% data from the natural sample for comparison.
%
%% Conversions and normalization

% Convert simulation time from y to Myr and extended volumes to radii.
simulationValuesXlTimeVsSize = [simulationValuesXlTimeVsSize(:,1)/10^6 ...
    (3*simulationValuesXlTimeVsSize(:,7)/(4*pi)).^(1/3)];

% Normalize the measured zoning-radius values to the radius of the first
% simulation crystal.
radiusNormFactor = max(simulationValuesTimeVsRadius(:,2))/...
    max(rockValuesGrtZoning(:,1));
rockValuesGrtZoning(:,1) = rockValuesGrtZoning(:,1)*radiusNormFactor;
rockValuesGrtZoningNorm = [rockValuesGrtZoning(:,1) rockValuesGrtZoning(:,2)];

%% Use Spline fits to derive a function for zoning vs time

% Determine the spline coefficients for the Measured zoning values as a
% function of radius.
rockZoningSplineCoef = ...
    spline(rockValuesGrtZoningNorm(:,1),rockValuesGrtZoningNorm(:,2));

% Use the spline coefficients to generate zoning values at each radius
% reported for the first crystal in the simulation, smooth the results, and
% associate time with the zoning values.
rockZoningFromRadius = ...
    ppval(rockZoningSplineCoef,simulationValuesTimeVsRadius(:,2));
rockZoningFromRadiusSmoothed = ...
    smooth(rockZoningFromRadius,smoothingWindowZoning);
simulationTimeAndZoning = ...
    [simulationValuesTimeVsRadius(:,1) rockZoningFromRadiusSmoothed];

% Check that time values are appropriate for the spline (below).  This has 
% benn problematic in past calculations.
if simulationTimeAndZoning(2,1) - simulationTimeAndZoning(1,1) == 0
    % If the first two rows are the same, something's wrong
    disp('---------------------------------------------');
    disp('   Error in calcZoningVsRadius.m');
    disp('   Time increment must increase in simulationTimeAndZoning');
    disp(['   There might be a problem with the ...radiiVsTime.txt file']);
    disp('---------------------------------------------');
    return
elseif simulationTimeAndZoning(end,1) - simulationTimeAndZoning((end-1),1) == 0
    % Remove duplicate entry in last row
    simulationTimeAndZoning = simulationTimeAndZoning(1:(end-1),:);
end

% Determine the spline coefficients for the Simulation zoning values as a
% function of time.
simulationTimeAndZoningSplineCoef = ...
    spline(simulationTimeAndZoning(:,1),simulationTimeAndZoning(:,2));

% Use the spline coefficients to generate zoning values at each nucleation 
% time reported from the simulation, and smooth the results.
simulationXlZoningFromTime = ...
    ppval(simulationTimeAndZoningSplineCoef,simulationValuesXlTimeVsSize(:,1));
simulationZoningFromTimeSmoothed = ...
    smooth(simulationXlZoningFromTime,smoothingWindowZoning);
simulationRadiusAndCoreConc = ...
    [simulationValuesXlTimeVsSize(:,2) simulationZoningFromTimeSmoothed];

end