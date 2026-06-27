%% main attitude determination script

clear;
close all

format longG

addpath(genpath('./SCT/SCTAcademic'));

%% mission parameters

% get orbital elements at some epoch
% use issorbit for initial model
% TODO: add polar orbit parameters
% TODO: update with austronauts mission parameters once available
% el    (1,6)   Elements vector [a,i,W,w,e,M]
% jD0   (1,1)   Julian date of epoch
[el, jD0] = ISSOrbit('fixed');

% get current pos + vel at orbital elements (ECI)
[r_sat, v_sat] = El2RV(el); % [km, km/s]

% initialise attitude quaternion to LVLH (ECI -> LVLH)
% note: SCT puts scalar at q1
q_sat = QLVLH(r_sat, v_sat);


% state vector
% x = [rx; ry; rz; vx; vy; vz; q1; q2; q3; q4]
x = [r_sat; v_sat; q_sat];


% visualise orbit path and locations
PltOrbit(el, jD0);
hold on
plot3(r_sat(1), r_sat(2), r_sat(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off

%% reference vectors (ECI)

% magnetic field reference vector
% use dipole for initial model
% TODO: improve with IGRF model
[b_ref, bDot_ref] = BDipole(r_sat, jD0, v_sat); % [T, T/s]

% convert magnetic reference vec to unit vector
uB_ref = b_ref./norm(b_ref);


% sun reference vector
% TODO: improve with SunV2 and compare with SunVectorECI
% u         (3,:)   Unit sun vector (vector TO the sun)
% r         (1,:)   Distance from origin to sun (km)
[uS_ref, rS_ref] = SunV1(jD0, r_sat);

% hold on
% arrow3(r_sat', (r_sat+uSun_ref)', 'y0.1')
% hold off

%% sensor models (body frame)

uB_meas = MagModel(x, jD0);

uS_meas = CSSModel(x, uS_ref);

%% TRIAD 

% TODO: investigate choice of first vector (see Wertz pg 425 and footnote)

A_TRIAD = TRIAD([uS_meas, uB_meas], [uS_ref, uB_ref]);
q_TRIAD = Mat2Q(A_TRIAD);

%% Attitude error

q_err = QMult(q_sat, quatinv(q_TRIAD')')

% renormalise
q_err = q_err/norm(q_err);

% theta
theta_err = 2*acos(q_err(1))
% in deg
rad2deg(theta_err)