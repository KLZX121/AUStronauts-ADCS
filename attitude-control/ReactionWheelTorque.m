
% function T = ReactionWheelTorque(q_d, q_c)

%% Program Summary

% INPUTS
% q_desired             (4, 1) desired direction quaternion input (default 
% from RHSCubesat is ISS facing earth and is calculated using QLVLH on the 
% ISS position and velocity)

% q_current             (4, 1) current direction quaternion input (default 
% from RHSCubeSat is [1,0,0,0] and is stored in x(7:10))

% OUTPUT
% T                     (3, 1) torque vector (x,y,z)

%% CubeSat Physical Data

d = RHSCubeSat; % initialising cubesat data

%%

model = '3U';

[area,nFace,rFace] = CubeSatFaces( model, 1 ); % get face data 
% (areas, normals, distance from centre)

% surface data
d.surfData.area = area;
d.surfData.nFace = nFace;
d.surfData.rFace = rFace;
d.surfData.att.type = 'eci';
d.atm = [];

d.mass = 6; % kg
d.inertia = InertiaCubeSat(model, d.mass); % assuming uniform distribution
% TODO: update to our cubesat data from structure team (mass, inertia,
% areas)

d.kWheels = 14:16; % indices of reaction wheels


%% ISS Model of orbit for current direction [WILL CHANGE DEPENDING
% ON OUR ORBIT]

x = d.x0; % default
[el, jD0] = ISSOrbit;
[r,v] = El2RV(el);
x(1:3) = r; % position vector
x(4:6) = v; % velocity vector
d.jD0 = jD0; % starting Julian date

%% Find desired rotation angle

% to get default desired and current quaternions
% if (nargin == 0)
    q_current = x(7:10);
    q_desired = QLVLH(x(1:3),x(4:6));
% else % input
%     q_current = q_c;
%     q_desired = q_d;
% end 

% Find the delta quaternion (change in quaternion needed)
delta_q     = QMult( QPose(q_current), q_desired);
[angle, u]  = Q2AU( delta_q );

% the angle is the amount of rotation needed around the unit vector u
% calculate a single torque value based on the desired rotation

time = 60; 
% time wanted to complete torque maneuver, can be changed and will result 
% in very different torques.

angular_acceleration = 2*angle / time^2;

single_torque = d.inertia * angular_acceleration; %(1, 1)

d.TRW = single_torque*u;