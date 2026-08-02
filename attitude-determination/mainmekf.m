clear;
close all;

if isempty(which('QProd'))
    addpath(genpath('../common'));
end

%% specs

% gyro variance
sigv = sqrt(10)*1e-7;
sigu = sqrt(10)*1e-10;
sigs = 0;
sigU = 0;
sigL = 0;

% gyro bias, scale factors, misalignments
bTrue0 = deg2rad([0.1; 0.1; 0.1])./3600;
sTrue = 1e-6*[1500; 1000; 1500];
kUTrue = 1e-6*[1000; 1500; 2000];
kLTrue = 1e-6*[500; 1000; 1500];


% magnetometer 1-sigma noise (cubemag compact)
d.mag.sigma = repmat((120/3)*1e-9, 3, 1);
d.mag.bias = zeros(3, 1);
d.mag.D = zeros(3, 1);
d.mag.O = eye(3);

qTrue0 = (sqrt(2)/2).*[1; 1; 0; 0];

% truth functions

wTrueFn = @(t) deg2rad(10).*[sin(0.01*t); sin(0.0085*t); cos(0.0085*t)];
bTrueFn = @(bOld, dt) bOld + sigu*sqrt(dt)*randn(3, 1);

%% initialise

mekfNom = MEKFNominal;
x0 = zeros(6, 1);
q0 = qTrue0;
P0 = blkdiag( ...
    deg2rad(6/3600)^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3) ...
);
STrue = MEKF.SMatrix(sTrue, kUTrue, kLTrue);
Q = blkdiag(sigv^2*eye(3), sigu^2*eye(3));
mekfNom.Initialise(x0, q0, P0, STrue, Q, d.mag);


mekfCal = MEKFGyroCal;
x0C = zeros(15, 1);
P0C = blkdiag( ...
    P0, ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3) ...
);
QC = blkdiag(sigv^2*eye(3), sigu^2*eye(3), sigs^2*eye(3), sigU^2*eye(3), sigL^2*eye(3));
mekfCal.Initialise(x0C, q0, P0C, STrue, QC, d.mag);


qTrue = qTrue0;
bTrue = bTrue0;
wTrue = wTrueFn(0);

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

xErrList = zeros(mekfNom.nStates, nSim+1);
xSigma3List = zeros(mekfNom.nStates, nSim+1);
xSigma3List(:, 1) = mekfNom.CalcSigma3();

xErrListC = zeros(mekfCal.nStates, nSim+1);
xSigma3ListC = zeros(mekfCal.nStates, nSim+1);
xSigma3ListC(:, 1) = mekfCal.CalcSigma3();

tic
for i = 1:nSim
    wTrue = wTrueFn(t(i));

    % system state
    r = rOrb(:, i);
    jD = jDOrb(i);

    wGyro = MEKF.GyroMeasurement(wTrue, STrue, bTrue, sigv/sqrt(dt));
    bRef = IGRFECI(r*1e3, jD);

    mekfNom.Step(dt, bRef, qTrue, wGyro);
    mekfCal.Step(dt, bRef, qTrue, wGyro);

    % propagate
    bTrue = bTrueFn(bTrue, dt);
    qTrue = MEKF.PropQDisc(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);


    xErrList(:, i+1) = mekfNom.CalcError(qTrue, bTrue);
    xSigma3List(:, i+1) = mekfNom.CalcSigma3();
    xErrListC(:, i+1) = mekfCal.CalcError(qTrue, [bTrue; sTrue; kUTrue; kLTrue]);
    xSigma3ListC(:, i+1) = mekfCal.CalcSigma3();
end
toc

mekfNom.Plot(t, xErrList, xSigma3List);
mekfCal.Plot(t, xErrListC, xSigma3ListC);