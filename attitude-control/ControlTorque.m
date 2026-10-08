function T_c = ControlTorque(q_d, x, d)

%% Program Summary

%   INPUTS
%   q_d             (4, 1)  desired direction quaternion input (default 
% from RHSCubesat is ISS facing earth and is calculated using QLVLH on the 
% ISS position and velocity)
%   time            (1, 1)  time allowed for procedure, e.g. 60s
%   x               (:, 1)  state vector
%   d                         Input parameters (struct)
%       .TExt       (3, 1)  Total external applied torque (N m)
%       .TRW        (3, 1)  Total reaction wheel torque (N m)
%       .ISat       (3, 3)  Moment of Inertia tensor of satellite (kg m^2)
%       .IRW        (1, 1)  Moment of Inertia of Reaction Wheels (kg m^2)
%       .iWSat      (1, 3)  Indices of satellite ang vel in state vector
%       .iWRW       (1, 3)  Indices of rw ang vel in state vector
%
% OUTPUT
% T                     (3, 1) control torque vector (x,y,z)

%% CubeSat Data

inertiaSat = InertiaCubeSat('3U', 6); % assuming uniform distribution
% TODO: update to our cubesat data from structure team (mass, inertia,
% areas)

%% Find delta quaternion

% input
q_current = x(d.iQ);
q_desired = q_d;
 
% Find the error quaternion (change in quaternion needed)
delta_q = QProd( q_current, QConj(q_desired));

%% Control torque from control

k_p = 1.62e-4; % proportional gain (CHANGE)
k_d = 9.75e-4; % derivative gain (CHANGE)
w = x(d.iWSat); % angular velocity of satellite

% Calculate the control torque vector
%T_c = -k_p * sign(delta_q(1)) * delta_q(2:4) - k_d * w;
t_c = -k_p * sign(delta_q(1)) * delta_q(2:4) - k_d * (1 +- delta_q(2:4).' * delta_q(2:4)) * w;

% limit to 586 rad/s
for j = 1:3
    if x(d.iWRW(j)) > 586
        %x(d.iWRW(j)) = 586;
        T_c(j) = 0;
    elseif x(d.iWRW(j)) < -586
        %x(d.iWRW(j)) = -586;
        T_c(j) = 0;
    else
        T_c(j) = -t_c(j);
    end
end

end