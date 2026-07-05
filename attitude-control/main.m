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
% attitude quaternion
[el, jD0] = ISSOrbit('fixed');
[rSat, vSat] = El2RV(el); % [km, km/s]
q0 = QLVLH(rSat, vSat);

%% plant

% TODO: wrap values somewhere else (another file?)
% TODO: define body frame (sct has z-axis longitudinal)
d.I = InertiaCubeSat('3U', 6);
d.iQ = 1:4;
d.iOmega = 5:7;
d.T = [1 1 1]';

% state vector
% TODO: properly define a state vector
x0 = [
    q0;
    omega0;
];

% ode state vector functions
% TODO: better define our data structure and integrate it less messily
xDotFn = @(x, t, d) [
    QKinematics(x(d.iQ), x(d.iOmega)); 
    EulerDynamics(x(d.iOmega), d);
];


%% propogate with integrator (RK4)
% TODO: compare RK4 with ode45 or other integrators
% TODO: write own RK4 to better suit the structure of our code
h = 0.001;
t0 = 0;
tf = 1;

t = t0:h:tf;
N = (tf-t0)/h;

xList = zeros(N+1, size(x0, 1));
xList(1, :) = x0';

x = x0;
for i = 1:N
    x = RK4(xDotFn, x, h, t(i), d);
    xList(i+1, :) = x';
end

%out = RK4Convergence(x0, xDotFn, d, 15);

%% plot
figure('Name', 'State Variables');
tiledlayout(2, 1)

nexttile
plot(t, xList(:, d.iQ))
legend("qs", 'q2', 'q3', 'q4')
grid on
ylabel('q')

nexttile
plot(t, xList(:, d.iOmega))
legend('wx', 'wy', 'wz')
grid on
ylabel('\omega (rad/s)')
xlabel('t (s)')

%% analyse attitude with euler axis/angle
% convert attitude quaternion into reference with initial quaternion
% TODO: QMult's syntax is contradictory - write our own function
qBL = QMult(xList(:, d.iQ)', QPose(q0));
qs = qBL(1, :)';
qv = qBL(2:4, :)';
eulAngle = 2*acos(qs);
eulAxis = qv./vecnorm(qv, 2, 2);

figure('Name', 'Attitude (LVLH -> body)');
tiledlayout(3, 1)

nexttile
plot(t, qBL)
legend('qs', 'q2', 'q3', 'q4')
ylabel('quaternions')
grid on

nexttile
plot(t, rad2deg(eulAngle))
ylabel('euler angle (deg)')
grid on

nexttile
plot(t, eulAxis)
ylabel('euler axis')
grid on
legend('x', 'y', 'z')

AnimQ(qBL)