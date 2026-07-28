clear;
close all;

%% specs

% variance
sigv = sqrt(10)*1e-7;
sigu = sqrt(10)*1e-10;

% cubemag compact 1-sigma noise
sigmaMag = (120/3)*1e-9;
RMag = sigmaMag^2*eye(3);

% bias, scale factors, misalignments
bTrue0 = deg2rad([0.1; 0.1; 0.1])./3600;
sTrue = 1e-6*[1500; 1000; 1500];
kUTrue = 1e-6*[1000; 1500; 2000];
kLTrue = 1e-6*[500; 1000; 1500];

qTrue0 = (sqrt(2)/2).*[1; 1; 0; 0];

P0 = blkdiag( ...
    deg2rad(6/3600)^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3) ...
);

%% truth functions

wTrueFn = @(t) deg2rad(10).*[sin(0.01*t); sin(0.0085*t); cos(0.0085*t)];
bTrueFn = @(bOld, dt) bOld + sigu*sqrt(dt)*randn(3, 1);

qTrue = qTrue0;
bTrue = bTrue0;
wTrue = wTrueFn(0);

%% initialise

x0 = zeros(6, 1);
q0 = qTrue0;

% calibration matrix
SFn = @(s, kU, kL) diag(s) + [0 kU(1) kU(2); kL(1) 0 kU(3); kL(2) kL(3) 0];
STrue = SFn(sTrue, kUTrue, kLTrue);

Q = blkdiag(sigv^2*eye(3), sigu^2*eye(3));
R = RMag;



mekf = MissionMEKF;
mekf.Initialise(x0, q0, P0, STrue, Q, R);

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

xErrList = zeros(6, nSim+1);
xSigma3List = zeros(6, nSim+1);
xSigma3List(:, 1) = mekf.CalcSigma3();

for i = 1:nSim
    wTrue = wTrueFn(t(i));

    extData.jD = jDOrb(i);
    extData.r = rOrb(:, i);
    extData.qTrue = qTrue;
    extData.sigmaMag = sigmaMag;
    extData.wGyro = MEKF.GyroMeasurement(wTrue, STrue, bTrue, sigv/sqrt(dt));


    mekf.Step(dt, extData);


    % propagate
    bTrue = bTrueFn(bTrue, dt);
    qTrue = MEKF.PropQDisc(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);


    xErrList(:, i+1) = mekf.CalcError(qTrue, bTrue);
    xSigma3List(:, i+1) = mekf.CalcSigma3();
end

mekf.Plot(t, xErrList, xSigma3List);