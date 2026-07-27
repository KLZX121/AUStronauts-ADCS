function T_c = ControlTorque(q_d, q_c, time)

%% Program Summary

% INPUTS
% q_d                   (4, 1) desired direction quaternion input (default 
% from RHSCubesat is ISS facing earth and is calculated using QLVLH on the 
% ISS position and velocity)

% q_c                   (4, 1) current direction quaternion input (default 
% from RHSCubeSat is [1,0,0,0] and is stored in x(7:10))

% time                  (1, 1) time allowed for procedure, e.g. 60s

% OUTPUT
% T                     (3, 1) control torque vector (x,y,z)

%% CubeSat Data

inertiaSat = InertiaCubeSat('3U', 6); % assuming uniform distribution
% TODO: update to our cubesat data from structure team (mass, inertia,
% areas)

%% Find desired rotation angle

% input
q_current = q_c;
q_desired = q_d;
t = time;
 
% Find the delta quaternion (change in quaternion needed)
delta_q     = QProd( QConj(q_current), q_desired);
[angle, u]  = QAngleUnit( delta_q );

% the angle is the amount of rotation needed around the unit vector u
% calculate a single torque value based on the desired rotation

% time wanted to complete torque maneuver, can be changed and will result 
% in very different torques.

angular_acceleration = 2*angle / t^2;

single_torque = inertiaSat * angular_acceleration; %(1, 1)

T_c = single_torque*u; % torque * unit vector

end