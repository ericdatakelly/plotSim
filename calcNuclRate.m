function [temperature,rate] = ...
    calcNuclRate(values,CiAl,k1,k2,GtoC,Teq)

% calcNuclRate.m
% This function calculates the nucleation rate at each time step
% using the following equation...
% dNv/dt = k1 * exp((-k2 * (CiAl)^GtoC * (Teq+273.15))/(delC^GtoC * (T+273.15)))



temperature = values(:,3);
delC = values(:,4)-CiAl;
rate = k1.* exp(-(k2*(CiAl^GtoC)*(Teq+273.15))./...
    ((delC.^GtoC).*(temperature+273.15)));
end
