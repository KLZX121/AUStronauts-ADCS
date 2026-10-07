%% main attitude control script

clear;
close all;

format longG
set(0, 'DefaultLegendLocation', 'eastoutside')
set(0, 'DefaultLineLineWidth', 1.4)
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultTextFontSize', 12)

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
if isempty(which('PropState'))
    addpath(genpath('../common'));
end

%% initial states

% orbital elements to find initial pos and vel
% TODO: switch to using meters and rewrite relevant functions
[el, jD0] = ISSOrbit('fixed');
[r0, v0] = El2RV(el); % [km, km/s]

% attitude quaternion
q0 = GetLVLHQ(r0, v0);

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

nStates = 16;


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
d.TRW = [0; 0; 0;]; % N m
% max rw speed
d.maxWRW = 5600*2*pi/60; % rad/s

% indices of states in state vetor
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;
d.iWSat = 11:13;
d.iWRW = 14:16;

%% plant

% ode state vector functions
% TODO: replace orbit propagator (FOrbCart) with mission team values
xDotFn = @(x, d) [
    FOrbCart(x);
    QKinematics(x, d); 
    EulerDynamics(x, d);
    RWDynamics(EulerDynamics(x, d), d);
];

%% propagate with integrator (RK4)
h = 1;
t0 = 0;
tf = 90*60;

t = t0:h:tf;
nSim = length(t);


qLVLHList = zeros(4, nSim);
qLVLHList(:, 1) = GetLVLHQ(r0, v0);
xList = zeros(nStates, nSim);
xList(:, 1) = x0;
distList = zeros(12, nSim); % [Tg; Ta; Ts; Tm]

x = x0;

tic
for i = 2:nSim
    % compute disturbances
    % TODO: rewrite our own disturbance functions since these use sct's q
    TGrav = GravityGradientFromR(QConj(x(d.iQ)), d.ISat, x(d.iR), 3.98600436e5);

xList = zeros(length(x0), length(t));
xList(:, 1) = x0;

x = x0;
for i = 2:length(t)
    if mod(i,3) == 0 || i == 2
        T_c = ControlTorque([1; 0; 0; 0], x, d);
        
    end

    x = PropState(xDotFn, x, d, h, t(i-1):h:t(i));
    xList(:, i) = x;
end
    jD = jD0 + t(i)/86400;
    s = SunV1(jD, x(d.iR));
    B = BDipole(x(d.iR), jD);
    [Tg, Ta, Ts, Tm] = DisturbanceTorques(x, d, s, B);

    d.TExt = Tg + Ta + Ts + Tm;


    x = PropState(xDotFn, x, d, h);

    % renormalise quaternion
    x(d.iQ) = x(d.iQ) ./ norm(x(d.iQ));


    xList(:, i) = x;
    qLVLHList(:, i) = GetLVLHQ(x(d.iR), x(d.iV));
    distList(:, i) = [Tg; Ta; Ts; Tm];
end
toc

%d.TExt = zeros(3, 1);
%out = RK4Convergence(x0, xDotFn, d, 8, 6, struct("t0", t0, "tf", tf, "h0", 10));

%% plot

% convert time to minutes
t = t./60;

figure('Name', 'State Variables');
tl = tiledlayout(4, 2);
tl.Title.String = 'State Variables (ECI)';
tl.Title.FontWeight = 'bold';

nexttile
plot(t, xList(d.iR, :))
legend('r_x', 'r_y', 'r_z')
grid on
ylabel('r (km)')
ylim('padded')

nexttile
plot(t, xList(d.iV, :))
legend('v_x', 'v_y', 'v_z')
grid on
ylabel('v (km/s)')
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

xlabel('t (min)')

%% attitude change plots

% euler axis/angle
% convert attitude quaternion into reference with initial quaternion
qBL = zeros(4, nSim);

for i = 1:nSim
    qBL(:, i) = QProd(xList(d.iQ, i), QConj(qLVLHList(:, i)));
    % enforce sign continuity
    if ( i >  1 && (dot(qBL(:, i), qBL(:, i-1)) < 0))
        qBL(:, i) = -qBL(:, i);
    end
end

% TODO: put this in a function
qs = qBL(1, :);
qv = qBL(2:4, :);
eulAngle = 2*acos(min(qs, 1));
eulAxis = qv./vecnorm(qv, 2, 1);

figure('Name', 'Attitude (LVLH -> body)');
tiledlayout(3, 1)

nexttile
plot(t, qBL)
legend('q_s', 'q_x', 'q_y', 'q_z')
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

%{
eulAngles = zeros(3, nSim);
for i = 1:nSim
    % TODO: rewrite Q2Eul
    eulAngles(:, i) = Q2Eul(qBL(:, i));
end

nexttile
plot(t, rad2deg(eulAngles))
grid on
legend('x', 'y', 'z')
ylabel('euler angles (deg)')
xlabel('t (min)')
%}

%AnimQ(qBL);

%% disturbance plots

figure('Name', 'Disturbance Torques (Body)')
tiledlayout(4, 1)

nexttile
plot(t, distList(1:3, :))
title('Gravity Gradient')
ylabel('T_g (N m)')

legend('x', 'y', 'z')

nexttile
plot(t, distList(4:6, :))
title('Aerodynamic Drag')
ylabel('T_a (N m)')

nexttile
plot(t, distList(7:9, :))
title('SRP')
ylabel('T_s (N m)')

nexttile
plot(t, distList(10:12, :))
title('Magnetic Field')
ylabel('T_m (N m)')
xlabel('t (min)')


