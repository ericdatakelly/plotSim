function [timeD,tempD,D] = calcD(Dinf,phi,tortuosity,Q,R,values)

% calcDeff.m
% This function calculates the diffusivity at the characteristic
% temperature, D(Tc).  A diffusivity of D(Tc) will produce the same
% cumulative diffusivity for the duration of the crystallization event as
% the integral of the D(T) taken across the duration of the event.


D1 = Dinf.*phi.*tortuosity.*(exp(-Q./(R.*(values(:,3)+273.15)))).*0.0001; % Convert to m^2/s by multiplying by 0.0001
D2 = D1.*values(2,2);
D2sum = sum(D2);
D3 = max(values(:,2)).*D1;
D4 = D2sum - D3;
Dindex = find(D4<0,1,'First');

D = D1(Dindex);
timeD = values(Dindex,2)/1000000;  % convert to millions of years
tempD = values(Dindex,3);
