%% main attitude control script

clear;
close all;

format longG

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
if isempty(which('PropState'))
    addpath(genpath('../common'));
end

%% initial states

% angular velocities
wSat0 = [0; 0; 0;]; % (rad/s)
wRW0 = [0; 0; 0;];  % (rad/s)

% attitude quaternion
[el, jD0] = ISSOrbit('fixed');
[rSat, vSat] = El2RV(el); % [km, km/s]
q0 = QLVLH(rSat, vSat);

%% plant

% data struct
% TODO: wrap values somewhere else (another file?)
% TODO: define body frame (sct has z-axis longitudinal)
d.ISat = InertiaCubeSat('3U', 6);   % satellite moi
d.TExt = [0 0 0]';            % external torques

d.IRW = (0.6e-3)/(5600*2*pi/60);    % rw moi
d.TRW = [0.2e-3; 0.2e-3; 0.2e-3;];                 % rw torques
% indices of states in state vetor
d.iQ = 1:4;
d.iWSat = 5:7;
d.iWRW = 8:10;

% state vector
% TODO: properly define a state vector
x0 = [
    q0;
    wSat0;
    wRW0;
];

% ode state vector functions
% TODO: better define our data structure and integrate it less messily
xDotFn = @(x, t, d) [
    QKinematics(x, d); 
    EulerDynamics(x, d);
    RWDynamics(EulerDynamics(x, d), d);
];


%% propagate with integrator (RK4)
% TODO: compare RK4 with ode45 or other integrators
% TODO: write own RK4 to better suit the structure of our code
h = 0.1;
t0 = 0;
tf = 60;
t = t0:h:tf;

[x, xList] = PropState(xDotFn, x0, d, h, t);

%out = RK4Convergence(x0, xDotFn, d, 15);

%% plot
figure('Name', 'State Variables');
tiledlayout(2, 1)

nexttile
plot(t, xList(d.iQ, :))
legend("qs", 'q2', 'q3', 'q4')
grid on
ylabel('q')

nexttile
plot(t, xList(d.iWSat, :))
legend('wx', 'wy', 'wz')
grid on
ylabel('\omega (rad/s)')
xlabel('t (s)')

%% attitude change plots
% euler axis/angle
% convert attitude quaternion into reference with initial quaternion
% TODO: QMult's syntax is contradictory - write our own function
qBL = QMult(xList(d.iQ, :), QPose(q0));
qs = qBL(1, :);
qv = qBL(2:4, :);
eulAngle = 2*acos(qs);
eulAxis = qv./vecnorm(qv, 2, 1);

figure('Name', 'Attitude (LVLH -> body)');
tiledlayout(4, 1)

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
ylim('padded')

% euler angles (3-2-1)

n = size(qBL, 2);

eulAngles = zeros(3, n);
for i = 1:n
    eulAngles(:, i) = Q2Eul(qBL(:, i));
end

nexttile
plot(t, rad2deg(eulAngles))
grid on
legend('x', 'y', 'z')
ylabel('euler angles (deg)')
xlabel('t (s)')


AnimQ(qBL);