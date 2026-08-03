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
q_desired = q_d; % [0.295347131961402;0.909299901960774;0.229703139008898;-0.182154406660284]; 
 
% Find the delta quaternion (change in quaternion needed)
delta_q     = QProd( q_desired, QConj(q_current));

%% Control torque from control

k_p = 1e-5; % proportional gain (CHANGE)
k_d = 0.5e-5; % derivative gain (CHANGE)
w = x(d.iWSat); % angular velocity of satellite

% Calculate the control torque vector
T_c = -k_p * sign(delta_q(1)) * delta_q(2:4) - k_d * w;


%% Comparison Torque

[angle, u]  = QAngleUnit( delta_q );
% the angle is the amount of rotation needed around the unit vector u
% calculate a single torque value based on the desired rotation

t = 60; % placeholder parameter

angular_acceleration = 2*angle / t^2;

single_torque = inertiaSat * angular_acceleration; %(1, 1)

T_c2 = single_torque*u; % torque * unit vector

end