%% main attitude control script

clear;
close all;

format longG
set(0, 'DefaultLegendLocation', 'eastoutside')
set(0, 'DefaultLineLineWidth', 1.4)

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
if isempty(which('PropState'))
    addpath(genpath('../common'));
end

%% initial states

% orbital elements to find initial pos and vel
[el, jD0] = ISSOrbit('fixed');
[r0, v0] = El2RV(el); % [km, km/s]

% attitude quaternion
q0 = QLVLH(r0, v0);

% angular velocities
wSat0 = [0; 0; 0;]; % (rad/s)
wRW0 = [0; 0; 0;];  % (rad/s)


% state vector
x0 = [
    r0;
    v0;
    q0;
    wSat0;
    wRW0;
];


% data struct
% TODO: define in separate file?
% TODO: define body frame (sct has z-axis longitudinal)
% satellite moi
d.ISat = InertiaCubeSat('3U', 6);
% external torques
d.TExt = [0 0 0]';
% rw moi
d.IRW = (0.6e-3)/(5600*2*pi/60);
% rw torques
d.TRW = [0.2e-3; 0.2e-3; 0.2e-3;];
% indices of states in state vetor
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;
d.iWSat = 11:13;
d.iWRW = 14:16;

%% plant

% ode state vector functions
% TODO: replace orbit propagator (FOrbCart) with mission team values
xDotFn = @(x, t, d) [
    FOrbCart(x);
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
tl = tiledlayout(4, 2);
tl.Title.String = 'State Variables (ECI)';
tl.Title.FontWeight = 'bold';

nexttile
plot(t, xList(d.iR, :))
legend('r_x', 'r_y', 'r_z')
grid on
ylabel('r (m)')
ylim('padded')

nexttile
plot(t, xList(d.iV, :))
legend('v_x', 'v_y', 'v_z')
grid on
ylabel('v (m/s)')
ylim('padded')

nexttile([2 1])
plot(t, xList(d.iQ, :))
legend('q_s', 'q_x', 'q_y', 'q_z')
grid on
ylabel('q')
ylim('padded')

nexttile([2 1])
plot(t, xList(d.iWSat, :))
legend('\omega_x', '\omega_y', '\omega_z')
grid on
ylabel('\omega_s_a_t (rad/s)')
ylim('padded')

nexttile
plot(t, xList(d.iWRW, :))
legend('\omega_x', '\omega_y', '\omega_z')
grid on
ylabel('\omega_r_w (rad/s)')
ylim('padded')

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


%AnimQ(qBL);