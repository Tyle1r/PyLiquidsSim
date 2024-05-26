function [M] = prandtlMeyerMforNu(nu,gamma)

% This function accepts the input argument nu in radians
% and solves for M that satisfies the prandtl-meyer relation
% Call it like this
% [M] = prandtlMeyerMforNu(nu,gamma)

const = sqrt((gamma+1)/(gamma-1));

% This is the symbolic prandtly meyer function which returns
% its result in radians
M = 1.5;
for iter = 1:100
    NuOfM = const*atan(sqrt(M^2-1)/const) -atan(sqrt(M^2-1));
    Residual = NuOfM - nu;
    dNudM = M/((M^2 - 1)^(1/2)*(((M^2 - 1)*(gamma - 1))/(gamma + 1) + 1)) - 1/(M*(M^2 - 1)^(1/2));
    M = M - Residual/dNudM;
    
    if (abs(Residual) < 1.0e-5)
        break;
    end
end