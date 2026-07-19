%% main attitude determination script

clear;
close all

format longG

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
% [r, v] = (3, n), t = (1, n)
[rOrb, vOrb, tOrb] = RVFromKepler(el); % [km, km/s]

%% initial state

% initial position and velocity
r0 = rOrb(:, 1);
v0 = vOrb(:, 1);

% initialise attitude quaternion to LVLH (ECI -> LVLH)
q0 = GetLVLHQ(r0, v0);

% initial state vector
x = [r0; v0; q0];

%% simulation loop
nSim = length(tOrb);

qList = zeros(4, nSim);
thetaErrList = zeros(nSim, 1);
qErrList = zeros(4, nSim);

for i = 1:nSim
    %%% update state
    r = rOrb(:, i);
    v = vOrb(:, i);
    jD = jD0 + tOrb(i)/86400;

    x = [r; v; q0;];


    %%% reference vectors (ECI)
    
    % magnetic field reference vector
    % use dipole for initial model
    % TODO: improve with IGRF model
    [bRef, bDotRef] = BDipole(r, jD, v); % [T, T/s]
    % convert to unit vector
    uBRef = bRef./norm(bRef);
    
    % sun reference vector
    % TODO: improve with SunV2 and compare with SunVectorECI
    [uSRef, rSRef] = SunV1(jD, r);
    

    %%% sensor models (body frame)
    
    uBMeas = MagModel(x, jD);
    uSMeas = CSSModel(x, uSRef);
    

    %%% TRIAD 
    
    % TODO: investigate choice of first vector (see Wertz pg 425 and footnote)
    
    ATRIAD = TRIAD([uBMeas, uSMeas], [uBRef, uSRef]);
    qTRIAD = DCMToQ(ATRIAD);

    %%% Attitude error
    [thetaErr, qErr] = QAttErr(q0, qTRIAD);


    %%% save
    qList(:, i) = qTRIAD;
    thetaErrList(i) = thetaErr;
    qErrList(:, i) = qErr;
end

%% plots
% visualise orbit path and positions
% TODO: visualise ref/body frames
% TODO: visualise reference/body vectors
PltOrbit(el, jD0);
hold on
plot3(rOrb(1, 1), rOrb(2, 1), rOrb(3, 1), 'or', 'MarkerSize', 10, 'LineWidth', 2)
hold off

figure('Name', 'Estimated Attitude')
plot(tOrb, qList, 'LineWidth', 1)
title('Estimated Quaternion')
ylabel('q')
xlabel('t (s)')
legend('q_s', 'q_x', 'q_y', 'q_z')
grid on

figure('Name', 'Attitude Error')
tiledlayout(2, 1)

nexttile
plot(tOrb, qErrList, 'LineWidth', 1)
title('Error Quaternion')
ylabel('q error')
xlabel('t (s)')
grid on
ylim('padded')
legend('q_s', 'q_x', 'q_y', 'q_z')

nexttile
plot(tOrb, rad2deg(thetaErrList), 'x-', 'MarkerSize', 8, 'LineWidth', 0.8)
title('Angular Error')
ylabel('\theta error (deg)')
xlabel('t (s)')
grid on



%Anim2Q([repmat(q0, 1, nSim); qList;])