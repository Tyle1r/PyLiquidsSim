function [f] = ffactor(Re,eps_D)

% Call this function like this
% f = ffactor(Re,epsilon/D)
%
% f is the pipe friction factor
% Re is the Reynolds number of the pipe flow
% epsilon/D is the roughness height divided by pipe diameter

% Just do a fixed point iteration
f = 0.02;

for i=1:100
    fnew = 1/(-2*log10(eps_D/3.7 +2.51/(Re*sqrt(f))))^2;
    if (abs(f-fnew) < 1.0e-8)
        break;
    end
    f = fnew;
end
if i > 99
    disp('Trouble converging to friction factor')
end