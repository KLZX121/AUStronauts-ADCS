clear;
close all;

if isempty(which('QProd'))
    addpath(genpath('../common'));
end
if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end

%%
% simulation timeframe and timestep
t0 = 0;
tf = 90*60;
dt = 1;
nSim = (tf-t0)/dt;
t = t0:dt:tf;

% orbit propagation
[el, jD0] = ISSOrbit('fixed');
[rOrb, vOrb] = RVFromKepler(el, t);
rOrb = rOrb.*1e3;
vOrb = vOrb.*1e3;
jDOrb = jD0:(dt/86400):(jD0+(tf/86400));

%% specs

% gyro
calG.bias = deg2rad([0.1; 0.1; 0.1])./3600;
calG.s = 1e-6*[1500; 1000; 1500];
calG.kU = 1e-6*[1000; 1500; 2000];
calG.kL = 1e-6*[500; 1000; 1500];
calG.sigG = repmat(sqrt(10)*1e-7, 3, 1);
calG.sigB = repmat(sqrt(10)*1e-10, 3, 1);

gyro = GyroModel(calG);

% magnetometer
calM.bias = zeros(3, 1);
calM.D = zeros(3, 3);
calM.O = eye(3);
calM.sigma = repmat((120/3)*1e-9, 3, 1);

mag = MagModel(calM);

% css
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
calCSS.sigmaDark = 1e-3;
calCSS.yLims = [0 2.4];
calCSS.T0 = 60;
calCSS.alpha = 4.31e-3;

css = CSSModel(calCSS);

% truth functions

T = 60;
qTrue = (sqrt(2)/2).*[1; 1; 0; 0];
wTrueFn = @(t) deg2rad(10).*[sin(0.01*t); sin(0.0085*t); cos(0.0085*t)];

%% initialise

bRef = IGRFECI(rOrb(:, 1), jD0);
uBRef = bRef./norm(bRef);

uSRef = SunV1(jD0, rOrb(:, 1));

bMag = mag.Measurement(qTrue, bRef, true);
uBMag = bMag/norm(bMag);

ySMeas = css.Measurement(qTrue, uSRef, T, true);
uSMeas = css.CalcSunVec(ySMeas, T);

ATRIAD = TRIAD([uBMag, uSMeas], [uBRef, uSRef]);
qTRIAD = DCMToQ(ATRIAD);

q0 = qTRIAD;




mekfNom = MEKFNominal;
x0 = zeros(6, 1);
P0 = blkdiag( ...
    deg2rad(0.2)^2*eye(3), ...
    (0.2*pi/(3600*180))^2.*eye(3) ...
);
Q = diag([gyro.sigG.^2; gyro.sigB.^2]);
mekfNom.Initialise(x0, q0, P0, Q, gyro, mag, css);


mekfCal = MEKFGyroCal;
x0C = zeros(15, 1);
P0C = blkdiag( ...
    P0, ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3), ...
    (0.002/3)^2.*eye(3) ...
);
QC = blkdiag(Q, zeros(9));
mekfCal.Initialise(x0C, q0, P0C, QC, gyro, mag, css);

%% simulation

% data save
xErrList = zeros(mekfNom.nStates, nSim+1);
xSigma3List = zeros(mekfNom.nStates, nSim+1);
xSigma3List(:, 1) = mekfNom.CalcSigma3();

xErrListC = zeros(mekfCal.nStates, nSim+1);
xSigma3ListC = zeros(mekfCal.nStates, nSim+1);
xSigma3ListC(:, 1) = mekfCal.CalcSigma3();

wTrue = wTrueFn(0);
tic
for i = 1:nSim
    % system state
    r = rOrb(:, i);
    jD = jDOrb(i);

    wGyro = gyro.Measurement(wTrue, dt);
    bRef = IGRFECI(r, jD);

    mekfNom.Step(dt, qTrue, bRef, uSRef, T, wGyro);
    mekfCal.Step(dt, qTrue, bRef, uSRef, T, wGyro);

    % propagate
    gyro.PropagateBias(dt);
    qTrue = MEKF.PropQDisc(qTrue, wTrue, dt);
    qTrue = qTrue/norm(qTrue);
    wTrue = wTrueFn(t(i+1));

    xErrList(:, i+1) = mekfNom.CalcError(qTrue, gyro.bias);
    xSigma3List(:, i+1) = mekfNom.CalcSigma3();
    xErrListC(:, i+1) = mekfCal.CalcError(qTrue, [gyro.bias; gyro.s; gyro.kU; gyro.kL]);
    xSigma3ListC(:, i+1) = mekfCal.CalcSigma3();
end
toc

mekfNom.Plot(t, xErrList, xSigma3List);
mekfCal.Plot(t, xErrListC, xSigma3ListC);