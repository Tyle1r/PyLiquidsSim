function [value,isterminal,direction] = myEventsFcn(t,y)

% When this value crosses zero there is an event
value = y(1)-500e3;
% I don't care if it is increasing or decreasing
direction = 0;
% Force the simulation to stop
isterminal = 1;
end