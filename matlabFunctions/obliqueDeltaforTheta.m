function [delta] = obliqueDeltaforTheta(m1,theta,gamma)

% This function accepts the input arguments of M1, theta, and gamma
% M1 is the upstream mach number
% theta is the oblique shock angle relative to the incoming flow
% gamma is the specific heat ratio

% It returns the turning angle delta 
% call it like this
% [delta] = obliqueDeltaforTheta(m1,theta,gamma)


tanDelta = (2*cot(theta)*(m1^2*sin(theta)^2-1))/(m1^2*(gamma+cos(2*theta))+2);
delta = atan(tanDelta);