function [index1,index2] = calcFirstLastNucl(values)

% calcFirstLastNucl.m
% This function calculates the index of first nucleation and last nucleation,
% at 95% Al consumed.  The indices are used to determine the time and temperature
% of the crystallization range.
% **Note: The accuracy of the indexing depends on the reporting interval!**


% Find the first crystal
	if ((min(values(:,10)) > 0) && (min(values(:,10)) == 1))
        index1 = find(values(:,10) == 1, 1 );
	else
        index1 = find(values(:,10) > 0, 1 );
	end


% Find the last crystal (the nucleation at 95% Al consumed)
	productAl95 = 0.95*max(values(:,6)); % the amount of Al in the porphyroblasts when 95% of the Al has been consumed
	index2 = find(values(:,6) > productAl95, 1 );  % find the line with a value near 95% Al and use this for the index of the end of crystallization
