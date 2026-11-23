function T_c = ControlTorque(q_d, x, d)

%% Program Summary

%   INPUTS
%   q_d             (4, 1)  desired direction quaternion input (default 
% from RHSCubesat is ISS facing earth and is calculated using QLVLH on the 
% ISS position and velocity)
%   x               (:, 1)  state vector
%   d                         Input parameters (struct)
%       .iWSat      (1, 3)  Indices of satellite ang vel in state vector
%       .iWRW       (1, 3)  Indices of rw ang vel in state vector
%
%   OUTPUT
%   T_c             (3, 1) control torque vector (x,y,z)

%% Find delta quaternion

% input
q_current = x(d.iQ);
q_desired = q_d;
 
% Find the error quaternion (change in quaternion needed)
delta_q = QProd( q_current, QConj(q_desired));
w = x(d.iWSat); % current satellite angular velocities


%temporary fix, change input to include a mode (slew or track)
    % code to find amount of turn
    deg_q = 2*acos(abs(delta_q(1)/norm(delta_q(1:4))))*180/pi;
if deg_q < 5
    operating_mode = "small_delq";
else
    % Big delta q (large angle flip, move to LVLH, slew)
    operating_mode = "slew";
end


if operating_mode == "small_delq"

    k_p = 1.62e-4; % proportional gain (CHANGE)
    k_d = 9.75e-4; % derivative gain (CHANGE)
    % Calculate the control torque vector
    T_c = -k_p * sign(delta_q(1)) * delta_q(2:4) - k_d ...
        * (1 + delta_q(2:4).' * delta_q(2:4)) * w;

elseif operating_mode == "slew"
    
    % Calculate the control torque vector

    %step 1: find angle between quaternions: deg_q
    %step 2: choose interpolation coefficient: T
    T = 0.82/deg_q;

    %step 3: interpolate to interpolation coefficient using slerp
    q_int = sin((1-T)*deg_q*pi/180) / sin(deg_q*pi/180) * -q_current ...
        + sin(T*deg_q*pi/180) / sin(deg_q*pi/180) * q_desired;

    %step 4: find error quaternion to interpolated quaternion
    delta_q_int = QProd( q_current, QConj(q_int));

    %step 5: calculate the control torque vector
    k_p = 1.62e-4; % proportional gain (CHANGE)
    k_d = 9.75e-4; % derivative gain (CHANGE)
    T_c = -k_p * delta_q_int(2:4) - k_d * w;
    
    %deg_q_int = 2*acos(abs(delta_q_int(1)/norm(delta_q_int(1:4))))*180/pi %debug

end


% Limit to max RW speed

for j = 1:3
    if abs(x(d.iWRW(j))) > d.maxWRW
        T_c(j) = 0;
    end
end

%T_c %debug
end