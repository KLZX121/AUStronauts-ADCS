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
if (q0(1) < 0)
    q0 = q0.*-1;
end

% initial state vector
x = [r0; v0; q0];

% temperature
T = 60;

d.jD = jD0;
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;

%% sensor models


calM.bias = zeros(3, 1);
calM.D = zeros(3, 3);
calM.O = eye(3);
calM.sigma = (120/3)*1e-9;

mag = MagModel(calM);


calCSS.nSensors = 6;
calCSS.uSensors = [eye(3) -eye(3)];
calCSS.fov = deg2rad(120);
calCSS.LUT = [
    0 2.0737;
    10 2.0442;
    20 1.9592;
    30 1.8058;
    40 1.5854;
    50 1.3486;
    60 1.0999;
    ];
calCSS.LUT(:, 1) = deg2rad(calCSS.LUT(:, 1));
calCSS.sigmaTheta = deg2rad(5);
calCSS.sigmaEdge = deg2rad(10);
calCSS.sigmaDark = 1e-3;
calCSS.yLims = [0 2.4];
calCSS.T0 = 60;
calCSS.alpha = 4.31e-3;

css = CSSModel(calCSS);

%% simulation loop
nSim = length(tOrb);

qList = zeros(8, nSim);
thetaErrList = zeros(1, nSim);
qErrList = zeros(4, nSim);
bList = zeros(6, nSim);
ySList = zeros(6, nSim);
sList = zeros(6, nSim);
yLitList = zeros(6, nSim);

q = q0;
a = deg2rad(10).*randn(3, 1);
b = deg2rad(5).*randn(3, 1);
c = deg2rad(5)*randn(3, 1);
wFn = @(t) a.*sin(b.*t) + c;
w = wFn(0);

tic
for i = 1:nSim
    %%% update state
    r = rOrb(:, i);
    v = vOrb(:, i);
    jD = jD0 + tOrb(i)/86400;
    d.jD = jD;

    x = [r; v; q;];


    %%% reference vectors (ECI)
    
    % magnetic field reference vector
    bRef = IGRFECI(r, jD); % T
    uBRef = bRef./norm(bRef);
    
    % sun reference vector
    % TODO: improve with SunV2 and compare with SunVectorECI
    uSRef = SunV1(jD, r);
    

    %%% sensor models (body frame)
    
    bMag = mag.Measurement(x(d.iQ), bRef, true);
    uBMag = bMag/norm(bMag);

    [ySMeas, yLit] = css.Measurement(x(d.iQ), uSRef, T, true);
    uSMeas = css.CalcSunVec(ySMeas, T);
    

    %%% TRIAD 
    
    % TODO: investigate choice of first vector (see Wertz pg 425 and footnote)
    
    [ATRIAD, PTRIAD] = TRIAD([uBMag, uSMeas], [uBRef, uSRef], mag.sigma);
    qTRIAD = DCMToQ(ATRIAD);
    if (qTRIAD(1) < 0)
        qTRIAD = qTRIAD.*-1;
    end

    %%% Attitude error
    [thetaErr, qErr] = QAttErr(q, qTRIAD);


    %%% save
    qList(1:4, i) = x(d.iQ);
    qList(5:8, i) = qTRIAD;

    thetaErrList(i) = thetaErr;
    qErrList(:, i) = qErr;

    bList(1:3, i) = QToDCM(x(d.iQ))*bRef;
    bList(4:6, i) = bMag;

    ySList(:, i) = ySMeas;

    sList(1:3, i) = QToDCM(x(d.iQ))*uSRef;
    sList(4:6, i) = uSMeas;

    yLitList(:, i) = yLit;


    %%% propagate
    q = MEKF.PropQDisc(q, w, tOrb(2) - tOrb(1));
    q = q/norm(q);
    if (q(1) < 0)
        q = q.*-1;
    end
    if (i < nSim)
        w = wFn(tOrb(i+1));
    end
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
figure('Name', 'Magnetometer');
tl = tiledlayout(3, 1);
tl.Title.String = "Magnetic Fields (Body)";
tl.Title.FontWeight = "bold";

bList = bList.*1e9;

nexttile
plot(tOrb./60, bList(1, :), '-')
hold on
plot(tOrb./60, bList(4, :), 'x')
hold off
grid on
xticklabels({})
ylabel('b_x (nT)')
ylim('padded')

legend('IGRF', 'Mag')

nexttile
plot(tOrb./60, bList(2, :), '-')
hold on
plot(tOrb./60, bList(5, :), 'x')
hold off
grid on
xticklabels({})
ylabel('b_y (nT)')
ylim('padded')

nexttile
plot(tOrb./60, bList(3, :), '-')
hold on
plot(tOrb./60, bList(6, :), 'x')
hold off
grid on
ylabel('b_z (nT)')
ylim('padded')

xlabel('t (min)')

% plot sun measurements
figure('Name', 'Sun Sensors');
tl = tiledlayout(4, 1);
tl.Title.FontWeight = "bold";

nexttile
colororder(lines(6))
col = lines(6);
box on
hold on
for i = 1:nSim
    for j = 1:css.nSensors
        litStatus = yLitList(j, i);

        if litStatus == 1 % in fov
            spec = 'o';
        elseif litStatus == 2 % outside fov
            spec = 'x';
        elseif litStatus == 0 % unlit
            spec = '.';
        end

        plot(tOrb(i)./60, ySList(j, i), spec, 'Color', col(j, :))
    end
end
hold off
title('Sun Sensor Measurements')
grid on
xticklabels({})
ylabel('y_s (V)')
ylim('padded')
hold on
legP{1} = plot(nan);
legP{2} = plot(nan);
legP{3} = plot(nan);
legP{4} = plot(nan);
legP{5} = plot(nan);
legP{6} = plot(nan);

legP{7} = plot(nan, 'o', 'Color', 'k');
legP{8} = plot(nan, 'x', 'Color', 'k');
legP{9} = plot(nan, '.', 'Color', 'k');
hold off
legend([legP{:}], {'+x', '+y', '+z', '-x', '-y', '-z', 'lit in FOV', 'lit outside FOV', 'unlit'})

nexttile
plot(tOrb./60, sList(1, :), '-')
title("Sun Vectors (Body)")
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
figure('Name', 'Estimated Attitude')
colororder(lines(4))

plot(tOrb./60, qList(1:4, :), '-', 'DisplayName', 'q_t_r_u_e')
hold on
plot(tOrb./60, qList(5:8, :), 'x','DisplayName', 'q_T_R_I_A_D')
hold off
title('Estimated Quaternion (TRIAD)')
ylabel('q')
xlabel('t (min)')
legend()
grid on

% attitude error plots
figure('Name', 'Attitude Error')
tiledlayout(2, 1)

nexttile
plot(tOrb./60, qErrList)
title('Error Quaternion')
ylabel('q error')
xlabel('t (min)')
grid on
ylim('padded')
legend('q_s', 'q_x', 'q_y', 'q_z')

nexttile
plot(tOrb./60, rad2deg(thetaErrList), '-', 'MarkerSize', 8)
title('Angular Error')
ylabel('\theta error (deg)')
xlabel('t (min)')
grid on

hold on
indices = find(sum(yLitList == 1) == 3);
for i = 1:length(indices)
    xline(tOrb(indices(i))./60)
end
hold off

%Anim2Q([repmat(q0, 1, nSim); qList;])







