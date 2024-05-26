function [thetas] = obliqueThetasforDelta(m1,delta,gamma)

% This function accepts the input arguments of M1, delta, and gamma
% and solves for the strong and weak shong angles (theta) that satisfies the oblique shock relation
% that give you a specific turning angle
% Call it like this
%  [thetas] = obliqueThetasforDelta(m1,delta,gamma)

tandelta = tan(delta);
% Newton Iteration to solve equation for weak shock
theta = asin(1/m1);
for iter = 1:100
    Eq = (2*cot(theta)*(m1^2*sin(theta)^2-1))/(m1^2*(gamma+cos(2*theta))+2);
    Err = Eq -tandelta;
    dErrdTheta = (4*m1^2*sin(2*theta)*cot(theta)*(m1^2*sin(theta)^2 - 1))/((cos(2*theta) + gamma)*m1^2 + 2)^2 - (2*(cot(theta)^2 + 1)*(m1^2*sin(theta)^2 - 1))/((cos(2*theta) + gamma)*m1^2 + 2) + (4*m1^2*cos(theta)*cot(theta)*sin(theta))/((cos(2*theta) + gamma)*m1^2 + 2);
    theta = theta - Err/dErrdTheta;
    theta = abs(theta);
    if (abs(Err) < 1.0e-5)
        break;
    end
end
if iter > 99 || theta < asin(1/m1) || theta > pi/2
    disp(['did not converge to weak shock solution, M = ' num2str(m1) ' \delta = ' num2str(delta*180/pi)]);
    thetas(1) = nan;
else
    thetas(1) = theta;
end

%Newton Iteration to solve equation for strong shock
theta = pi/2;
for iter = 1:100
    Eq = (2*cot(theta)*(m1^2*sin(theta)^2-1))/(m1^2*(gamma+cos(2*theta))+2);
    Err = Eq -tandelta;
    dErrdTheta = (4*m1^2*sin(2*theta)*cot(theta)*(m1^2*sin(theta)^2 - 1))/((cos(2*theta) + gamma)*m1^2 + 2)^2 - (2*(cot(theta)^2 + 1)*(m1^2*sin(theta)^2 - 1))/((cos(2*theta) + gamma)*m1^2 + 2) + (4*m1^2*cos(theta)*cot(theta)*sin(theta))/((cos(2*theta) + gamma)*m1^2 + 2);
    theta = theta - Err/dErrdTheta;
    theta = abs(theta);
    if (abs(Err) < 1.0e-5)
        break;
    end
end
if iter > 99 || theta < asin(1/m1) || theta > pi/2
    disp(['did not converge to strong shock solution, M = ' num2str(m1) ' \delta = ' num2str(delta*180/pi)])
    thetas(2) = nan;
else
    thetas(2) = theta;
end

if thetas(1) > thetas(2)
    disp('strong and weak are inverted?')
    thetas
end