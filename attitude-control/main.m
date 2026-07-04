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
x0 = [omega0];

%% plant

% ode state vector functions
% TODO: better define our data structure and integrate it less messily
xDotFn = [@(x, t, d) EulerDynamics(x(1:3), d)];

% TODO: wrap values somewhere else (another file?)
% TODO: define body frame (sct has z-axis longitudinal)
d.I = InertiaCubeSat('3U', 6);
d.T = [1 1 1]';


% propogate with integrator (RK4)
% TODO: compare RK4 with ode45 or other integrators
% TODO: write own RK4 to better suit the structure of our code
h = 0.001;
t0 = 0;
tf = 1;

t = t0:h:tf;
N = (tf-t0)/h;

x = x0;
for i = 1:N
    x = RK4(xDotFn, x, h, t(i), d);
end

x

out = RK4Convergence(x0, xDotFn, d, 15);