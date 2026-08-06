clear;
close all

format longG
set(0, 'DefaultLegendLocation', 'eastoutside')
set(0, 'DefaultLineLineWidth', 1.4)
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultTextFontSize', 12)

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
if isempty(which('QProd'))
    addpath(genpath('../common'));
end



lit = zeros(3, 2);
tic
for i = 1:1000
    i
    lit = lit + run;
end
toc

% calc mean
lit(:, 2) = lit(:, 2) ./ lit(:, 1);

nRes = sum(lit, 1);
nRes = nRes(1);

% convert n to %
lit(:, 1) = lit(:, 1) ./ nRes;
% convert err to deg
lit(:, 2) = rad2deg(lit(:, 2));

lit


function lit = run

[el, jD0] = ISSOrbit('fixed');

% randomise orbit inclination
el(2) = el(2) + deg2rad(10)*randn;


[rOrb, vOrb, tOrb] = RVFromKepler(el);
rOrb = rOrb.*1e3;
vOrb = vOrb.*1e3;

% initial state

% initial position and velocity
r0 = rOrb(:, 1);
v0 = vOrb(:, 1);

% initialise attitude quaternion to LVLH (ECI -> LVLH)
q0 = GetLVLHQ(r0, v0);

% initial state vector
x = [r0; v0; q0];

% temperature
T = 60;

d.jD = jD0;
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;

% sensor models


calM.bias = zeros(3, 1);
calM.D = zeros(3, 3);
calM.O = eye(3);
calM.sigma = repmat((120/3)*1e-9, 3, 1);

mag = MagModel(calM);


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
calCSS.sigmaEdge = deg2rad(10);
calCSS.sigmaDark = 1e-3;
calCSS.yLims = [0 2.4];
calCSS.T0 = 60;
calCSS.alpha = 4.31e-3;

css = CSSModel(calCSS);

% simulation loop
nSim = length(tOrb);

qList = zeros(8, nSim);
thetaErrList = zeros(1, nSim);
qErrList = zeros(4, nSim);
bList = zeros(6, nSim);
ySList = zeros(6, nSim);
sList = zeros(6, nSim);
yLitList = zeros(6, nSim);

q = q0;
wFn = @(t) deg2rad(10).*randn(3, 1).*sin(0.1.*randn(3, 1).*t) + deg2rad(5)*randn(3, 1);
w = wFn(0);

for i = 1:nSim
    %%% update state
    r = rOrb(:, i);
    v = vOrb(:, i);
    jD = jD0 + tOrb(i)/86400;
    d.jD = jD;

    x = [r; v; q;];


    %%% reference vectors (ECI)
    
    % magnetic field reference vector
    bRef = IGRFECI(r, jD); % T
    uBRef = bRef./norm(bRef);
    
    % sun reference vector
    % TODO: improve with SunV2 and compare with SunVectorECI
    uSRef = SunV1(jD, r);
    

    %%% sensor models (body frame)
    
    bMag = mag.Measurement(x(d.iQ), bRef, true);
    uBMag = bMag/norm(bMag);

    [ySMeas, yLit] = css.Measurement(x(d.iQ), uSRef, T, true);
    uSMeas = css.CalcSunVec(ySMeas, T);
    

    %%% TRIAD 
    
    % TODO: investigate choice of first vector (see Wertz pg 425 and footnote)
    
    ATRIAD = TRIAD([uBMag, uSMeas], [uBRef, uSRef]);
    qTRIAD = DCMToQ(ATRIAD);
    if (qTRIAD(1) < 0)
        qTRIAD = qTRIAD.*-1;
    end

    %%% Attitude error
    [thetaErr, qErr] = QAttErr(q, qTRIAD);


    %%% save
    qList(1:4, i) = x(d.iQ);
    qList(5:8, i) = qTRIAD;

    thetaErrList(i) = thetaErr;
    qErrList(:, i) = qErr;

    bList(1:3, i) = QToDCM(x(d.iQ))*bRef;
    bList(4:6, i) = bMag;

    ySList(:, i) = ySMeas;

    sList(1:3, i) = QToDCM(x(d.iQ))*uSRef;
    sList(4:6, i) = uSMeas;

    yLitList(:, i) = yLit;

    %%% propagate
    q = MEKF.PropQDisc(q, w, tOrb(2) - tOrb(1));
    q = q/norm(q);
    if (q(1) < 0)
        q = q.*-1;
    end
    if (i < nSim)
        w = wFn(tOrb(i+1));
    end
end










% 3 lit fov
indices = find(sum(yLitList == 1) == 3);
n3Lit = length(indices);
sumlit3 = 0;
for i = 1:n3Lit
   sumlit3 = sumlit3 + thetaErrList(indices(i));
end

% 2 lit fov
indices = find(sum(yLitList == 1) == 2);
n2Lit = length(indices);
sumlit2 = 0;
for i = 1:n2Lit
    sumlit2 = sumlit2 + thetaErrList(indices(i));
end

% 1 lit fov
indices = find(sum(yLitList == 1) == 1);
n1Lit = length(indices);
sumlit1 = 0;
for i = 1:n1Lit
    sumlit1 = sumlit1 + thetaErrList(indices(i));
end

lit = [
    n3Lit sumlit3;
    n2Lit sumlit2;
    n1Lit sumlit1;
];
end