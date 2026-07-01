%% main attitude control script

clear;
close all;

format longG

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end

%% initial state
% angular velocities
omega0 = [0; 0; 0;]; % (rad/s)

% state vector
% TODO: properly define a state vector
x = [omega0];

%% plant

% ode state vector functions
% TODO: better define our data structure and integrate it less messily
xDotFn = [@(x, t, T, I) EulerDynamics(x(1:3), T, I)];

% TODO: wrap values somewhere else (another file?)
% TODO: define body frame (sct has z-axis longitudinal)
I = InertiaCubeSat('3U', 6);

T = [1; 1; 1;];

% propogate with integrator (RK4)
% TODO: compare RK4 with ode45 or other integrators
% TODO: write own RK4 to better suit the structure of our code
h = 5e-5;
t0 = 0;
tf = 1;

t = t0:h:tf;
N = (tf-t0)/h;


for i = 1:N
    x = RK4(xDotFn, x, h, t(i), T, I);
end

x
