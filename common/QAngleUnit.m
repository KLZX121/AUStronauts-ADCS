function [angle, u] = QAngleUnit(q)

% QAngleUnit calculates the angle and unit vector of a quaternion
%   [angle, u] = QAngleUnit(q)
%
%   Inputs
%   q       (4, 1) quaternion to convert
%
%   Outputs
%   angle   (1, 1) angle the satellite must turn
%   u       (3, 1) unit vector at which the angle turns around

angle = 2*acos( q(1));

if angle < -pi
    angle = angle + 2*pi;
elseif angle > pi
    angle = angle - 2*pi;
end

if angle ~= 0 % remove div 0 errors
    u = q(2:4)/sin(angle/2);
else
    u = [1,0,0];
end

end