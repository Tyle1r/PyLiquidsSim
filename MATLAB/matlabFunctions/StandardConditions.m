function [T,P,rho] = StandardConditions(altitude)
% Reads data listed in Table 2.9 of book  and interpolates at altitude
% altitude should be in meters
% The function should be in the current folder or MATLAB Search Path along
% with the file std.csv


persistent std_data;
if isempty(std_data)
    [std_data] = csvread('std.csv');
    std_data(:,3) = log(std_data(:,3));
    std_data(:,4) = log(std_data(:,4));
end

% Use interp1 to interpolate the values in the table

[data] = interp1(std_data(:,1),std_data(:,2:end),altitude);
T = data(:,1);
P = 1.0133e5*exp(data(:,2));
rho = exp(data(:,3));
end