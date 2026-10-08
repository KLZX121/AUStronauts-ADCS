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

% controller frequency
d.controlHZ = 1;

% indices of states in state vetor
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;
d.iWSat = 11:13;
d.iWRW = 14:16;

% satellite geometry
geometry = jsondecode(fileread('../geometry.json'));

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


results.qLVLH = zeros(4, nSim);
results.qLVLH(:, 1) = GetLVLHQ(r0, v0);
results.x = zeros(nStates, nSim);
results.x(:, 1) = x0;
results.TC = zeros(3, nSim);
results.dist(1:nSim) = struct( ...
    'g', zeros(3, 1), ...
    'a', zeros(3, 1), ...
    's', zeros(3, 1), ...
    'm', zeros(3, 1) ...
);

surfForce.a.f = zeros(3, geometry.n_surfaces);
surfForce.a.t = surfForce.a.f;
surfForce.s.f = surfForce.a.f;
surfForce.s.t = surfForce.a.f;
results.surfForces(1:nSim) = surfForce;

results.s = zeros(3, nSim);
results.rs = zeros(1, nSim);
[results.s(:, 1), results.rs(1)] = SunV1(jD0, r0);

x = x0;

tic
for i = 2:nSim
    % lvlh
    qLVLH = GetLVLHQ(x(d.iR), x(d.iV));

    % compute control torque
    if i == 2 || mod(i*h,d.controlHZ) == 0
        TC = ControlTorque(qLVLH, x, d);
        d.TRW = -TC;
    end

    % compute disturbances
    % TODO: get environmental densities
    % TODO: get srp
    jD = jD0 + t(i)/86400;
    [env.s, rs] = SunV1(jD, x(d.iR));
    env.B = BDipole(x(d.iR), jD);
    env.rho = 3.8e-12; % kg m^-3
    env.p = 4.5e-6; % N m^-2

    [distT, distSurf] = DisturbanceTorques(x, d, env, geometry);
    d.TExt = distT.g + distT.a + distT.s + distT.m;


    x = PropState(xDotFn, x, d, h);

    % renormalise quaternion
    x(d.iQ) = x(d.iQ) ./ norm(x(d.iQ));


    results.x(:, i) = x;
    results.TC(:, i) = d.TRW;
    results.qLVLH(:, i) = qLVLH;
    results.dist(i) = distT;
    results.surfForces(i) = distSurf;
    results.s(:, i) = env.s;
    results.rs(i) = rs;
end
toc


%d.TExt = zeros(3, 1);
%out = RK4Convergence(x0, xDotFn, d, 8, 6, struct("t0", t0, "tf", tf, "h0", 10));


% write results to files
mkdir('../results');
writelines(jsonencode(results.x'), '../results/xdata.json')
for i = 1:nSim
    surfOut(i).a.f = results.surfForces(i).a.f';
    surfOut(i).a.t = results.surfForces(i).a.t';
    surfOut(i).s.f = results.surfForces(i).s.f';
    surfOut(i).s.t = results.surfForces(i).s.t';
end
writelines(jsonencode(surfOut), '../results/surfdata.json')
writelines(jsonencode((results.s .* results.rs)'), '../results/rs.json')

%% plot

% convert time to minutes
t = t./60;

figure('Name', 'State Variables');
tl = tiledlayout(4, 2);
tl.Title.String = 'State Variables (ECI)';
tl.Title.FontWeight = 'bold';

nexttile
plot(t, results.x(d.iR, :))
legend('r_x', 'r_y', 'r_z')
grid on
ylabel('r (km)')
ylim('padded')

nexttile
plot(t, results.x(d.iV, :))
legend('v_x', 'v_y', 'v_z')
grid on
ylabel('v (km/s)')
ylim('padded')

nexttile([2 1])
plot(t, results.x(d.iQ, :))
legend('q_s', 'q_x', 'q_y', 'q_z')
grid on
ylabel('q')
ylim('padded')

nexttile([2 1])
plot(t, results.x(d.iWSat, :))
legend('\omega_x', '\omega_y', '\omega_z')
grid on
ylabel('\omega_s_a_t (rad/s)')
ylim('padded')

nexttile
plot(t, results.x(d.iWRW, :))
legend('\omega_x', '\omega_y', '\omega_z')
grid on
ylabel('\omega_r_w (rad/s)')
ylim('padded')

nexttile
plot(t, results.TC)
title('Command Torques')
legend('T_x', 'T_y', 'T_z')
grid on
ylabel('T_r_w (N/m)')
ylim('padded')

xlabel('t (min)')

%% attitude change plots

% euler axis/angle
% convert attitude quaternion into reference with initial quaternion
qBL = zeros(4, nSim);

for i = 1:nSim
    qBL(:, i) = QProd(results.x(d.iQ, i), QConj(results.qLVLH(:, i)));
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
plot(t, [results.dist.g])
title('Gravity Gradient')
ylabel('T_g (N m)')
grid on

legend('x', 'y', 'z')

nexttile
plot(t, [results.dist.a])
title('Aerodynamic Drag')
ylabel('T_a (N m)')
grid on

nexttile
plot(t, [results.dist.s])
title('SRP')
ylabel('T_s (N m)')
grid on

nexttile
plot(t, [results.dist.m])
title('Magnetic Field')
ylabel('T_m (N m)')
xlabel('t (min)')
grid on