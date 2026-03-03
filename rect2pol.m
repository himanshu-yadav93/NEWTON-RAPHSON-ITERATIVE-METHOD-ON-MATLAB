% Polar to Rectangular Conversion
% X - Complex matrix or number, X = A + jB
% RHO - Magnitude
% THETA - Angle in radians

function [rho, theta] = rect2pol(x)
rho = abs(x);
theta = angle(x);