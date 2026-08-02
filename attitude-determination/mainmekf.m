clear;
close all;

if isempty(which('QProd'))
    addpath(genpath('../common'));
end

%% specs

% gyro variance
d.gyro.sigG = sqrt(10)*1e-7;
d.gyro.sigB = sqrt(10)*1e-10;

% gyro bias, scale factors, misalignments
d.gyro.bias = deg2rad([0.1; 0.1; 0.1])./3600;
d.gyro.s = 1e-6*[1500; 1000; 1500];
d.gyro.kU = 1e-6*[1000; 1500; 2000];
d.gyro.kL = 1e-6*[500; 1000; 1500];
d.gyro.S = MEKF.SMatrix(d.gyro.s, d.gyro.kU, d.gyro.kL);


% magnetometer 1-sigma noise (cubemag compact)
d.mag.sigma = repmat((120/3)*1e-9, 3, 1);
d.mag.bias = zeros(3, 1);
d.mag.D = zeros(3, 1);
d.mag.O = eye(3);
d.mag.M = eye(3) / (eye(3) + d.mag.D);

% truth functions

wTrueFn = @(t) deg2rad(10).*[sin(0.01*t); sin(0.0085*t); cos(0.0085*t)];
bTrueFn = @(bOld, dt) bOld + d.gyro.sigB*sqrt(dt)*randn(3, 1);

%% initialise
q0 = (sqrt(2)/2).*[1; 1; 0; 0];

mekfNom = MEKFNominal;
x0 = zeros(6, 1);
P0 = blkdiag( ...
    deg2rad(6/3600)^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3) ...
);
Q = blkdiag(d.gyro.sigG^2*eye(3), d.gyro.sigB^2*eye(3));
mekfNom.Initialise(x0, q0, P0, Q, d.gyro, d.mag);


mekfCal = MEKFGyroCal;
x0C = zeros(15, 1);
P0C = blkdiag( ...
    P0, ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3) ...
);
QC = blkdiag(d.gyro.sigG^2*eye(3), d.gyro.sigB^2*eye(3), zeros(9, 9));
mekfCal.Initialise(x0C, q0, P0C, QC, d.gyro, d.mag);

%% simulation

% simulation timeframe and timestep
t0 = 0;
tf = 90*60;
dt = 1;
nSim = (tf-t0)/dt;
t = t0:dt:tf;

% orbit propagation
[el, jD0] = ISSOrbit('fixed');
[rOrb, vOrb] = RVFromKepler(el, t);
jDOrb = jD0:(dt/86400):(jD0+(tf/86400));

% data save
xErrList = zeros(mekfNom.nStates, nSim+1);
xSigma3List = zeros(mekfNom.nStates, nSim+1);
xSigma3List(:, 1) = mekfNom.CalcSigma3();

xErrListC = zeros(mekfCal.nStates, nSim+1);
xSigma3ListC = zeros(mekfCal.nStates, nSim+1);
xSigma3ListC(:, 1) = mekfCal.CalcSigma3();

qTrue = q0;
wTrue = wTrueFn(0);
tic
for i = 1:nSim

    % system state
    r = rOrb(:, i);
    jD = jDOrb(i);

    wGyro = GyroModel(wTrue, d.gyro, dt);
    bRef = IGRFECI(r*1e3, jD);

    mekfNom.Step(dt, bRef, qTrue, wGyro);
    mekfCal.Step(dt, bRef, qTrue, wGyro);

    % propagate
    d.gyro.bias = bTrueFn(d.gyro.bias, dt);
    qTrue = MEKF.PropQDisc(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);
    wTrue = wTrueFn(t(i+1));

    xErrList(:, i+1) = mekfNom.CalcError(qTrue, d.gyro.bias);
    xSigma3List(:, i+1) = mekfNom.CalcSigma3();
    xErrListC(:, i+1) = mekfCal.CalcError(qTrue, [d.gyro.bias; d.gyro.s; d.gyro.kU; d.gyro.kL]);
    xSigma3ListC(:, i+1) = mekfCal.CalcSigma3();
end
toc

mekfNom.Plot(t, xErrList, xSigma3List);
mekfCal.Plot(t, xErrListC, xSigma3ListC);