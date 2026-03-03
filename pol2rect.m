% Polar to Rectangular Conversion
% RECT - Complex matrix or number, RECT = A + jB, A = Real, B = Imaginary
% RHO - Magnitude
% THETA - Angle in radians

function rect = pol2rect(rho,theta)
rect = rho.*cos(theta) + 1i*rho.*sin(theta);