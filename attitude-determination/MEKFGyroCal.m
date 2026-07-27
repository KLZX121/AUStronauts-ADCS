% Gyro Calibration MEKF based on Markley Sect 6.2.2

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
% scale factors and misalignment are assumed constant
sigs = 0;
sigU = 0;
sigL = 0;

% spectral density matrix
Q = blkdiag(sigv^2*eye(3), sigu^2*eye(3), sigs^2*eye(3), sigU^2*eye(3), sigL^2*eye(3));

%% dynamics

qTrue0 = (sqrt(2)/2).*[1; 0; 0; 1];
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
s0 = zeros(3, 1);
kU0 = zeros(3, 1);
kL0 = zeros(3, 1);

x0 = [
    dAng0;
    b0;
    s0;
    kU0;
    kL0;
];
xTrue0 = [
    dAng0;
    bTrue0;
    sTrue;
    kUTrue;
    kLTrue;
];

NStates = length(x0);

% initial covariance matrix
% how uncertain is the initial estimates?
% R = sigma^2*I_3

P0 = blkdiag( ...
    (6*pi/(3600*180))^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3) ...
);

%% state model

function [F, G] = StateMatrices(wGyro, wEst, SEst, bEst)
    U = [
        wGyro(2)-bEst(2), wGyro(3)-bEst(3), 0;
        0, 0, wGyro(3)-bEst(3);
        0, 0, 0
    ];
    L = [
        0, 0, 0;
        wGyro(1)-bEst(1), 0, 0;
        0, wGyro(1)-bEst(1), wGyro(2)-bEst(2);
    ];

    F = [
        -Skew(wEst), -(eye(3) - SEst), -diag(wGyro - bEst), -U, -L;
        zeros(12, 15)
    ];
    G = blkdiag(-(eye(3)-SEst), eye(3), eye(3), eye(3), eye(3));
end

%% observation model

% measurement sensitivity matrix
h = zeros(3, 1);
H = [eye(3) zeros(3, NStates-3)];

% measurement covariance
sigmaMeas = 6*pi/(180*3600);
RMeas = sigmaMeas^2*eye(3);

R = blkdiag(RMeas);

% simulate star tracker measurement
function y = SimST(qTrue, qEst, sigmaMeas)
    noise = sigmaMeas*randn(3, 1);
    dq = [1; 0.5*noise];
    dq = dq./norm(dq);
    qStar = QProd(dq, qTrue);

    qMeas = QProd(qStar, QConj(qEst));
    y = 2*qMeas(2:4)./qMeas(1);
end

%% discrete propagation functions

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
tf = 90*60;
dt = 1;

nSim = (tf-t0)/dt;
t = t0:dt:tf;

xList = zeros(NStates, nSim);
xErrList = zeros(NStates, nSim);
sigma3List = zeros(NStates, nSim);

for k = 1:nSim
    %%% calculate gain
    
    K = P*H'*inv(H*P*H'+R);
    

    %%% update
    
    % update covariance
    P = (eye(size(P, 1))-K*H)*P;
    
    % update estimated state with measurements
    xEst(1:3) = zeros(3, 1);
    y = SimST(qTrue, qEst, sigmaMeas);
    xEst = xEst + K*(y - h);
    
    % update estimated quaternion with error angle
    qEst = qEst + QKinematics([qEst; xEst(1:3)], struct('iQ', 1:4, 'iWSat', 5:7));
    qEst = qEst./norm(qEst);
    
    % get estimated states
    bEst = xEst(4:6);
    sEst = xEst(7:9);
    kUEst = xEst(10:12);
    kLEst = xEst(13:15);
    SEst = SFn(sEst, kUEst, kLEst);


    %%% propagate (discrete)
    
    wGyro = wGyroFn(wTrue, bTrue, dt);    

    % calculate state models
    wEst = (eye(3) - SEst)*(wGyro - bEst);
    [F, G] = StateMatrices(wGyro, wEst, SEst, bEst);

    P = PropPDiscrete(P, F, G, Q, dt, NStates);
    qEst = PropQDiscrete(qEst, wEst, dt);


    %%% propagate true values (only used to model measurements, won't be in real filter)
    
    bTrue = bTrueFn(bTrue, dt);
    qTrue = PropQDiscrete(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);

    wTrue = wTrueFn(t(k+1));

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
    xErrList(4:15, k) = xTrue(4:15) - xEst(4:15);
    sigma3List(:, k) = 3*sqrt(max(diag(P),0));
end

xList = [x0 xList];
xErrList = [xTrue0 xErrList];
sigma3List = [3*sqrt(max(diag(P0),0)) sigma3List];

%% plots

% error
figure;
tl = tiledlayout(3, 2);
tl.Title.String = "MEKF Error Results";

nexttile([1 2])
plot(t./60, 1e6*xErrList(1:3,:))
grid on
ylabel('\delta\vartheta (\murad)')
xlabel('t (min)')
ylim([-20 20])

lg = legend('x', 'y', 'z');
lg.Location = "eastoutside";

nexttile
plot(t./60, rad2deg(xErrList(4:6, :)).*3600)
grid on
ylabel('\beta (deg/h)')
xlabel('t (min)')
ylim([-0.02 0.02])

nexttile
plot(t./60, 1e6*xErrList(7:9, :))
grid on
ylabel('s (ppm)')
xlabel('t (min)')
ylim([-60 60])

nexttile
plot(t./60, 1e6*xErrList(10:12, :))
grid on
ylabel('k_U (ppm)')
xlabel('t (min)')
ylim([-60 60])

nexttile
plot(t./60, 1e6*xErrList(13:15, :))
grid on
ylabel('k_L (ppm)')
xlabel('t (min)')
ylim([-60 60])

figure;
plot(t./60, 1e6*xList(1:3,:))
grid on
ylabel('\delta\vartheta (\murad)')
xlabel('t (min)')
ylim([-20 20])