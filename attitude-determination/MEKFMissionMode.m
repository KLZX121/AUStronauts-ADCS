% mission mode MEKF based on Markley Sect 6.2.4
% parameters based on example 6.2
% assuming that true scale factors and misalignments are available

clear;
close all;

if isempty(which('QProd'))
    addpath(genpath('../common'));
end

%% gyro calibration parameters

% bias, scale factors, misalignments
bTrue0 = deg2rad([0.1; 0.1; 0.1])./3600;
sTrue = 1e-6*[1500; 1000; 1500];
kUTrue = 1e-6*[1000; 1500; 2000];
kLTrue = 1e-6*[500; 1000; 1500];

% calibration matrix
SFn = @(s, kU, kL) [s(1) kU(1) kU(2); kL(1) s(2) kU(3); kL(2) kL(3) s(3)];
STrue = SFn(sTrue, kUTrue, kLTrue);

% variance
sigv = sqrt(10)*1e-7;
sigu = sqrt(10)*1e-10;

% spectral density matrix
Q = blkdiag(sigv^2*eye(3), sigu^2*eye(3));

%% dynamics

qTrue0 = (sqrt(2)/2).*[1; 1; 0; 0];

wTrueFn = @(t) deg2rad(0.1).*[sin(0.01*t); sin(0.0085*t); cos(0.0085*t)];

% gyro bias random walk
bTrueFn = @(bOld, dt) bOld + sigu*sqrt(dt)*randn(3, 1);

% model gyro measurement
wGyroFn = @(wTrue, bTrue, dt) (eye(3) + STrue)*wTrue + bTrue + sigv/sqrt(dt)*randn(3, 1);

%% initial state

% initialise attitude to true initial quaternion
% initialise calibration estimates to 0
q0 = qTrue0;
dAng0 = zeros(3, 1);
b0 = zeros(3, 1);

x0 = [
    dAng0;
    b0;
];
xTrue0 = [
    dAng0;
    bTrue0;
];

NStates = length(x0);

% initial covariance matrix
% how uncertain is the initial estimates?
% R = sigma^2*I_3

P0 = blkdiag( ...
    deg2rad(6/3600)^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3) ...
);

%% state model

function [F, G] = StateMatrices(wEst, STrue)
    F = [
        -Skew(wEst), -(eye(3)-STrue);
        zeros(3, 6);
    ];
    G = blkdiag(-(eye(3)-STrue), eye(3));
end

%% observation model (magnetometer)

% cubemag compact 1-sigma noise
sigmaMag = (120/3)*1e-9;
RMag = sigmaMag^2*eye(3);

% simulate magnetometer measurement and measurement sensitivity matrix
function [y, h, H] = SimMag(r, jD, qTrue, qEst, sigmaMag)
    dMag = MeasMagnetometerEarth;
    dMag.jD = jD;
    dMag.kR = 1:3;
    dMag.kQ = 4:7;
    dMag.quantization = 1e-20;
    
    x = [r; qTrue;];

    bMeas = MagModel(x, dMag);
    bEst = QToDCM(qEst)*BDipole(r, jD);

    y = bMeas + sigmaMag*randn(3, 1);

    H = [Skew(bEst) zeros(3, 3)];
    h = bEst;
end

%% measurement covariance

R = blkdiag(RMag);

%% discrete propagation functions (rewrite with eq 6.93)

function PNew = PropPDiscrete(POld, F, G, Q, dt, NStates)
    A = [
        -F, G*Q*G';
        zeros(NStates, NStates), F';
    ].*dt;

    B = expm(A);

    Phi = B(NStates+1:end, NStates+1:end)';
    QK = Phi*B(1:NStates, NStates+1:end);

    PNew = Phi*POld*Phi' + QK;
end

function qNew = PropQDiscrete(qOld, wEst, dt)
    dw = 0.5*norm(wEst)*dt;

    psi = sin(dw)*wEst/norm(wEst);

    Theta = [
        cos(dw), -psi';
        psi, cos(dw)*eye(3)-Skew(psi);
    ];

    qNew = Theta*qOld;
end

%% filter

%%% initialise
P = P0;
xEst = x0;
qEst = q0;

% true values
qTrue = qTrue0;
bTrue = bTrue0;
wTrue = wTrueFn(0);

xTrue = xTrue0;


%%% simulation loop

t0 = 0;
tf = 270*60;
dt = 5;

nSim = (tf-t0)/dt;
t = t0:dt:tf;

%%% orbit propagation
[el, jD0] = ISSOrbit('fixed');
[rOrb, vOrb] = RVFromKepler(el, t);
jDOrb = jD0:(dt/86400):(jD0+(tf/86400));

xList = zeros(NStates, nSim);
xErrList = zeros(NStates, nSim);
sigma3List = zeros(NStates, nSim);

for k = 1:nSim
    wTrue = wTrueFn(t(k));

    %%% simulate measurements
    [y, h, H] = SimMag(rOrb(:, k), jDOrb(k), qTrue, qEst, sigmaMag);
    wGyro = wGyroFn(wTrue, bTrue, dt);


    %%% calculate gain
    
    K = P*H'*inv(H*P*H'+R);
    

    %%% update
    
    % update covariance
    P = (eye(size(P, 1))-K*H)*P;
    
    % explicit reset of attitude errors
    xEst(1:3) = zeros(3, 1);
    % update estimated state
    % note: biases are implicitly reset
    xEst = xEst + K*(y - h);
    
    % update estimated quaternion with error angle
    qEst = qEst + QKinematics([qEst; xEst(1:3)], struct('iQ', 1:4, 'iWSat', 5:7));
    qEst = qEst./norm(qEst);


    %%% propagate (discrete)  

    % calculate state models
    bEst = xEst(4:6);
    wEst = (eye(3) - STrue)*(wGyro - bEst);
    [F, G] = StateMatrices(wEst, STrue);

    P = PropPDiscrete(P, F, G, Q, dt, NStates);
    qEst = PropQDiscrete(qEst, wEst, dt);


    %%% propagate true values (only used to model measurements, won't be in real filter)
    
    bTrue = bTrueFn(bTrue, dt);
    qTrue = PropQDiscrete(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);

    %%% save data

    xList(:, k) = xEst;

    qErr = QProd(qTrue, QConj(qEst));
    qErr = qErr/norm(qErr);
    
    if qErr(1) < 0
        qErr = -qErr;
    end
    
    dAngErr = 2*qErr(2:4)/qErr(1);
    xErrList(1:3,k) = dAngErr;

    xTrue(4:6) = bTrue;
    xErrList(4:6, k) = xTrue(4:6) - xEst(4:6);
    sigma3List(:, k) = 3*sqrt(max(diag(P),0));
end

xList = [x0 xList];
xErrList = [xTrue0 xErrList];
sigma3List = [3*sqrt(max(diag(P0),0)) sigma3List];

%% plots

% error
figure;
tl = tiledlayout(2, 1);
tl.Title.String = "MEKF Results";

nexttile()
plot(t./60, rad2deg(xErrList(1:3,:)))
grid on
ylabel('\delta\vartheta (deg)')
xlabel('t (min)')
%ylim([-0.1 2e-1])

lg = legend('x', 'y', 'z');
lg.Location = "eastoutside";

nexttile
plot(t./60, rad2deg(xErrList(4:6, :)).*3600)
grid on
ylabel('\beta (deg/h)')
xlabel('t (min)')
%ylim([-0.2 0.2])