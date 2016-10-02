function [reactant,product,fluid] = calcAl(values)

% calcAl.m
% This function calculates the amount of aluminum in the reactants, products, and fluid


total = (values(:,5) + values(:,6) + values(:,7));
reactant = 100*(values(:,5)./total);
product = 100*(values(:,6)./total);
fluid = 100*(values(:,7)./total);