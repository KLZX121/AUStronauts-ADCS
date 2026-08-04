%% main attitude determination script

clear;
close all

format longG
set(0, 'DefaultLegendLocation', 'eastoutside')
set(0, 'DefaultLineLineWidth', 1.4)
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultTextFontSize', 12)

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
if isempty(which('QProd'))
    addpath(genpath('../common'));
end

%% orbital parameters

% get orbital elements at some epoch
% use issorbit for initial model
% TODO: add polar orbit parameters
% TODO: update with austronauts mission parameters once available
% el    (1,6)   Elements vector [a,i,W,w,e,M]
% jD0   (1,1)   Julian date of epoch
[el, jD0] = ISSOrbit('fixed');

% get pos + vel at orbital elements for one orbit (ECI)
[rOrb, vOrb, tOrb] = RVFromKepler(el); % [km, km/s]
% convert to m
rOrb = rOrb.*1e3;
vOrb = vOrb.*1e3;

%% initial state

% initial position and velocity
r0 = rOrb(:, 1);
v0 = vOrb(:, 1);

% initialise attitude quaternion to LVLH (ECI -> LVLH)
q0 = GetLVLHQ(r0, v0);

% initial state vector
x = [r0; v0; q0];

% temperature
T = 60;

d.jD = jD0;
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;

%% sensor models


magBias = zeros(3, 1);
magD = zeros(3, 3);
magO = eye(3);
magSigma = repmat((120/3)*1e-9, 3, 1);

mag = MagModel(magBias, magD, magO, magSigma);


nSunSensors = 6;
uSunSensors = [eye(3) -eye(3)];
fov = deg2rad(120);
LUT = [
    0 2.0737;
    10 2.0442;
    20 1.9592;
    30 1.8058;
    40 1.5854;
    50 1.3486;
    60 1.0999;
    ];
LUT(:, 1) = deg2rad(LUT(:, 1));
sigmaSun = 0;
sigmaTheta = deg2rad(5/3);
yLims = [0 2.4];
T0 = 60;
alpha = 4.31e-3;

css = CSSModel(nSunSensors, uSunSensors, fov, LUT, sigmaSun, sigmaTheta, yLims, T0, alpha);

%% simulation loop
nSim = length(tOrb);

qList = zeros(8, nSim);
thetaErrList = zeros(1, nSim);
qErrList = zeros(4, nSim);
bList = zeros(9, nSim);
sList = zeros(6, nSim);

tic
for i = 1:nSim
    %%% update state
    r = rOrb(:, i);
    v = vOrb(:, i);
    jD = jD0 + tOrb(i)/86400;
    d.jD = jD;

    x = [r; v; q0;];


    %%% reference vectors (ECI)
    
    % magnetic field reference vector
    bRef = IGRFECI(r, jD); % T
    uBRef = bRef./norm(bRef);
    
    % sun reference vector
    % TODO: improve with SunV2 and compare with SunVectorECI
    uSRef = SunV1(jD, r);
    

    %%% sensor models (body frame)
    
    bMag = mag.Measurement(bRef, x(d.iQ), true);
    uBMag = bMag/norm(bMag);

    ySMeas = css.Measurement(x(d.iQ), uSRef, T, true);
    uSMeas = css.CalcSunVec(ySMeas, T);
    

    %%% TRIAD 
    
    % TODO: investigate choice of first vector (see Wertz pg 425 and footnote)
    
    ATRIAD = TRIAD([uBMag, uSMeas], [uBRef, uSRef]);
    qTRIAD = DCMToQ(ATRIAD);

    %%% Attitude error
    [thetaErr, qErr] = QAttErr(q0, qTRIAD);


    %%% save
    qList(1:4, i) = x(d.iQ);
    qList(5:8, i) = qTRIAD;

    thetaErrList(i) = thetaErr;
    qErrList(:, i) = qErr;

    bList(1:3, i) = QToDCM(x(d.iQ))*BDipole(r.*1e-3, jD);
    bList(4:6, i) = QToDCM(x(d.iQ))*bRef;
    bList(7:9, i) = bMag;

    sList(1:3, i) = QToDCM(x(d.iQ))*uSRef;
    sList(4:6, i) = uSMeas;
end
toc

%% plots

% visualise orbit path and positions
% TODO: visualise ref/body frames
% TODO: visualise reference/body vectors
PltOrbit(el, jD0);
hold on
plot3(rOrb(1, 1)*1e-3, rOrb(2, 1)*1e-3, rOrb(3, 1)*1e-3, 'or', 'MarkerSize', 10, 'LineWidth', 2)
hold off

% plot magnetic measurement
figure('Name', 'Magnetometer Measurements');
tl = tiledlayout(3, 1);
tl.Title.String = "Magnetic Fields (Body)";
tl.Title.FontWeight = "bold";

nexttile
plot(tOrb./60, bList(1, :), ':')
hold on
plot(tOrb./60, bList(4, :), '-')
plot(tOrb./60, bList(7, :), 'x')
hold off
grid on
xticklabels({})
ylabel('b_x (T)')
ylim('padded')

legend('BDipole', 'IGRF', 'Mag')

nexttile
plot(tOrb./60, bList(2, :), ':')
hold on
plot(tOrb./60, bList(5, :), '-')
plot(tOrb./60, bList(8, :), 'x')
hold off
grid on
xticklabels({})
ylabel('b_y (T)')
ylim('padded')

nexttile
plot(tOrb./60, bList(3, :), ':')
hold on
plot(tOrb./60, bList(6, :), '-')
plot(tOrb./60, bList(9, :), 'x')
hold off
grid on
ylabel('b_z (T)')
ylim('padded')

xlabel('t (min)')

% plot sun measurements
figure('Name', 'Sun Sensor Measurements');
tl = tiledlayout(3, 1);
tl.Title.String = "Sun Vectors (Body)";
tl.Title.FontWeight = "bold";

nexttile
plot(tOrb./60, sList(1, :), '-')
hold on
plot(tOrb./60, sList(4, :), 'x')
hold off
grid on
xticklabels({})
ylabel('s_x')
ylim('padded')

legend('s_r_e_f', 's_c_s_s')

nexttile
plot(tOrb./60, sList(2, :), '-')
hold on
plot(tOrb./60, sList(5, :), 'x')
hold off
grid on
xticklabels({})
ylabel('s_y')
ylim('padded')

nexttile
plot(tOrb./60, sList(3, :), '-')
hold on
plot(tOrb./60, sList(6, :), 'x')
hold off
grid on
ylabel('s_z')
ylim('padded')

xlabel('t (min)')

% attitude plot
figure('Name', 'Estimated TRIAD Attitude')
colororder(lines(4))

plot(tOrb, qList(1:4, :), '-', 'DisplayName', 'q_t_r_u_e')
hold on
plot(tOrb, qList(5:8, :), 'x','DisplayName', 'q_T_R_I_A_D')
hold off
title('Estimated Quaternion (TRIAD)')
ylabel('q')
xlabel('t (s)')
legend()
grid on

% attitude error plots
figure('Name', 'Attitude Error')
tiledlayout(2, 1)

nexttile
plot(tOrb, qErrList)
title('Error Quaternion')
ylabel('q error')
xlabel('t (s)')
grid on
ylim('padded')
legend('q_s', 'q_x', 'q_y', 'q_z')

nexttile
plot(tOrb, rad2deg(thetaErrList), '-', 'MarkerSize', 8)
title('Angular Error')
ylabel('\theta error (deg)')
xlabel('t (s)')
grid on



%Anim2Q([repmat(q0, 1, nSim); qList;])