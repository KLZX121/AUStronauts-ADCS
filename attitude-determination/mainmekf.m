clear;
close all;

if isempty(which('QProd'))
    addpath(genpath('../common'));
end

%% specs

% variance
sigv = sqrt(10)*1e-7;
sigu = sqrt(10)*1e-10;
sigs = 0;
sigU = 0;
sigL = 0;

% magnetometer 1-sigma noise (cubemag compact)
sigmaMag = (120/3)*1e-9;
RMag = sigmaMag^2*eye(3);

% star tracker
sigmaST = 6*pi/(180*3600);
RST = sigmaST^2*eye(3);

% bias, scale factors, misalignments
bTrue0 = deg2rad([0.1; 0.1; 0.1])./3600;
sTrue = 1e-6*[1500; 1000; 1500];
kUTrue = 1e-6*[1000; 1500; 2000];
kLTrue = 1e-6*[500; 1000; 1500];

qTrue0 = (sqrt(2)/2).*[1; 1; 0; 0];

% truth functions

wTrueFn = @(t) deg2rad(10).*[sin(0.01*t); sin(0.0085*t); cos(0.0085*t)];
bTrueFn = @(bOld, dt) bOld + sigu*sqrt(dt)*randn(3, 1);

qTrue = qTrue0;
bTrue = bTrue0;
wTrue = wTrueFn(0);

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
R = RMag;
mekfNom.Initialise(x0, q0, P0, STrue, Q, R);


mekfCal = MEKFGyroCal;
x0C = zeros(15, 1);
q0C = qTrue0;
P0C = blkdiag( ...
    deg2rad(6/3600)^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3) ...
);
STrueC = STrue;
QC = blkdiag(sigv^2*eye(3), sigu^2*eye(3), sigs^2*eye(3), sigU^2*eye(3), sigL^2*eye(3));
RC = RST;
mekfCal.Initialise(x0C, q0C, P0C, STrueC, QC, RC);

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

extData.sigmaMag = sigmaMag;
extData.sigmaST = sigmaST;
for i = 1:nSim
    wTrue = wTrueFn(t(i));

    extData.jD = jDOrb(i);
    extData.r = rOrb(:, i);
    extData.qTrue = qTrue;
    extData.wGyro = MEKF.GyroMeasurement(wTrue, STrue, bTrue, sigv/sqrt(dt));


    mekfNom.Step(dt, extData);
    mekfCal.Step(dt, extData);

    % propagate
    bTrue = bTrueFn(bTrue, dt);
    qTrue = MEKF.PropQDisc(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);


    xErrList(:, i+1) = mekfNom.CalcError(qTrue, bTrue);
    xSigma3List(:, i+1) = mekfNom.CalcSigma3();
    xErrListC(:, i+1) = mekfCal.CalcError(qTrue, [bTrue; sTrue; kUTrue; kLTrue]);
    xSigma3ListC(:, i+1) = mekfCal.CalcSigma3();
end

mekfNom.Plot(t, xErrList, xSigma3List);
mekfCal.Plot(t, xErrListC, xSigma3ListC);