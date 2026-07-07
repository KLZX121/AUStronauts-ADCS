%% main attitude determination script

clear;
close all

format longG

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end

%% mission parameters

% get orbital elements at some epoch
% use issorbit for initial model
% TODO: add polar orbit parameters
% TODO: update with austronauts mission parameters once available
% el    (1,6)   Elements vector [a,i,W,w,e,M]
% jD0   (1,1)   Julian date of epoch
[el, jD0] = ISSOrbit('fixed');

% get current pos + vel at orbital elements (ECI)
[r0, v0] = El2RV(el); % [km, km/s]

% initialise attitude quaternion to LVLH (ECI -> LVLH)
% note: SCT puts scalar at q1
q0 = QLVLH(r0, v0);


% state vector
% x = [rx; ry; rz; vx; vy; vz; q1; q2; q3; q4]
x = [r0; v0; q0];


% visualise orbit path and positions
% TODO: visualise multiple positions
% TODO: visualise ref/body frames
% TODO: visualise reference/body vectors
PltOrbit(el, jD0);
hold on
plot3(r0(1), r0(2), r0(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off

%% reference vectors (ECI)

% magnetic field reference vector
% use dipole for initial model
% TODO: improve with IGRF model
[bRef, bDotRef] = BDipole(r0, jD0, v0); % [T, T/s]

% convert magnetic reference vec to unit vector
uBRef = bRef./norm(bRef);


% sun reference vector
% TODO: improve with SunV2 and compare with SunVectorECI
% u         (3,:)   Unit sun vector (vector TO the sun)
% r         (1,:)   Distance from origin to sun (km)
[uSRef, rSRef] = SunV1(jD0, r0);

% hold on
% arrow3(r_sat', (r_sat+uSun_ref)', 'y0.1')
% hold off

%% sensor models (body frame)

uBMeas = MagModel(x, jD0);
uSMeas = CSSModel(x, uSRef);

%% TRIAD 

% TODO: investigate choice of first vector (see Wertz pg 425 and footnote)

ATRIAD = TRIAD([uSMeas, uBMeas], [uSRef, uBRef]);
qTRIAD = Mat2Q(ATRIAD);

%% Attitude error

[thetaErr, qErr] = QAttErr(q0, qTRIAD)

Anim2Q([q0; qTRIAD])