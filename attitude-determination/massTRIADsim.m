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



n = 1000;

CSSTRIADSig(n)

function CSSTRIADSig(n)
% finds sig1 and sig2 for triad covariance
% only considers 3-fov css measurements (since thats what mekf will use)
    tSig1 = zeros(n, 1);
    tSig2 = zeros(n, 1);
    tsigmag = zeros(n, 1);

    tic
    for i = 1:n
        fprintf('%d\n',i)

        [~, ~, ~, sig1, sig2, sigmag] = RunOrbit;

        tSig1(i) = sig1;
        tSig2(i) = sig2;

        tsigmag(i) = sigmag;

    end
    toc

    tSig1 = tSig1(~isnan(tSig1));
    tSig2 = tSig2(~isnan(tSig2));
    tsigmag = tsigmag(~isnan(tsigmag));

    % [sigmag, sig1, sig2,]
    siglist = [tsigmag, tSig1, tSig2];

    rmssig = zeros(3, 1);
    meansig = zeros(3, 1);
    stdsig = zeros(3, 1);
    p95 = zeros(3, 1);
    for j = 1:3
        sig = siglist(:, j);

        rmssig(j) = sqrt(mean(sig.^2));

        meansig(j) = mean(sig);

        
        
        stdsig(j) = std(sig);
        p95(j) = prctile(sig,95);
    end

    results = table( ...
        rmssig, meansig,stdsig,p95, ...
        'VariableNames',{ ...
        'RMS', ...
        'Mean', ...
        'Std', ...
        'P95'},'RowNames', {'sigma_m', 'sigma_1', 'sigma_2'});

    disp(results)

end


function CSSFOVSense(n)
    lit = zeros(3,2);
    tErr = zeros(n,100);
    nLitAll = zeros(n,100);
    
    tic
    for i = 1:n
        fprintf('%d\n',i)
    
        [l,t,nLit] = RunOrbit;
    
        lit = lit+l;
        tErr(i,:) = t;
        nLitAll(i,:) = nLit;
    end
    toc
    
    
    
    sensorCounts = [3;2;1];
    
    percentage = zeros(3,1);
    meanError = zeros(3,1);
    stdError = zeros(3,1);
    medianError = zeros(3,1);
    p95Error = zeros(3,1);
    
    for j = 1:3
        mask = nLitAll == sensorCounts(j);
    
        % Select every error having this number of in-FOV sensors
        errors = rad2deg(tErr(mask));
    
        percentage(j) = 100*nnz(mask)/numel(tErr);
        meanError(j) = mean(errors);
        stdError(j) = std(errors);
        medianError(j) = median(errors);
        p95Error(j) = prctile(errors,95);
    end
    
    results = table( ...
        sensorCounts,percentage,meanError,stdError,medianError,p95Error, ...
        'VariableNames',{ ...
        'SensorsInFOV', ...
        'Percentage', ...
        'MeanErrorDeg', ...
        'StdErrorDeg', ...
        'MedianErrorDeg', ...
        'P95ErrorDeg'});
    
    disp(results)
end






function [lit,thetaErrList,nLitList, sig1, sig2, sigmag] = RunOrbit
[el, jD0] = ISSOrbit('fixed');

% randomise orbit inclination
el(2) = el(2) + deg2rad(10)*randn;

[rOrb, vOrb, tOrb] = RVFromKepler(el);
rOrb = rOrb.*1e3;
vOrb = vOrb.*1e3;

r0 = rOrb(:, 1);
v0 = vOrb(:, 1);

q0 = GetLVLHQ(r0, v0);
if (q0(1) < 0)
    q0 = q0.*-1;
end

T = 60;

d.jD = jD0;
d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;

calM.bias = zeros(3, 1);
calM.D = zeros(3, 3);
calM.O = eye(3);
calM.sigma = (120/3)*1e-9;

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

nSim = length(tOrb);

% [sig1; sig2; sigmag]
siglist = zeros(3, nSim);

qList = zeros(8, nSim);
thetaErrList = zeros(1, nSim);
qErrList = zeros(4, nSim);
bList = zeros(6, nSim);
ySList = zeros(6, nSim);
sList = zeros(6, nSim);
% lit
yLitList = zeros(6, nSim);

q = q0;
a = deg2rad(10).*randn(3, 1);
b = deg2rad(5).*randn(3, 1);
c = deg2rad(5)*randn(3, 1);
wFn = @(t) a.*sin(b.*t) + c;
w = wFn(0);

for i = 1:nSim
    r = rOrb(:, i);
    v = vOrb(:, i);
    jD = jD0 + tOrb(i)/86400;
    d.jD = jD;

    x = [r; v; q;];

    bRef = IGRFECI(r, jD);
    uBRef = bRef./norm(bRef);
    
    uSRef = SunV1(jD, r);
    
    bMag = mag.Measurement(x(d.iQ), bRef, true);
    uBMag = bMag/norm(bMag);

    [ySMeas, yLit] = css.Measurement(x(d.iQ), uSRef, T, true);
    uSMeas = css.CalcSunVec(ySMeas, T);
    
    [ATRIAD] = TRIAD([uBMag, uSMeas], [uBRef, uSRef]);
    qTRIAD = DCMToQ(ATRIAD);
    if (qTRIAD(1) < 0)
        qTRIAD = qTRIAD.*-1;
    end

    [thetaErr, qErr] = QAttErr(q, qTRIAD);


    % triad sig (only get 3-fov css measurements)
    if (sum(yLit(yLit == 1)) == 3)
        trueA = QToDCM(q);
        uBReal = trueA*(bRef./norm(bRef));
        uSReal = trueA*uSRef;
        a = cross(uBReal, uSReal);
        a = a./norm(a);
        dvb = uBMag - uBReal;
        dvs = uSMeas - uSReal;
        siglist(1, i) = abs(dot(a, dvb));
        siglist(2, i) = abs(dot(a, dvs));
        siglist(3, i) = mag.sigma(1)/norm(bMag);
    end


    qList(1:4, i) = x(d.iQ);
    qList(5:8, i) = qTRIAD;

    thetaErrList(i) = thetaErr;
    qErrList(:, i) = qErr;

    bList(1:3, i) = QToDCM(x(d.iQ))*bRef;
    bList(4:6, i) = bMag;

    ySList(:, i) = ySMeas;

    sList(1:3, i) = QToDCM(x(d.iQ))*uSRef;
    sList(4:6, i) = uSMeas;


    % lit
    yLitList(:, i) = yLit;


    q = MEKF.PropQDisc(q, w, tOrb(2) - tOrb(1));
    q = q/norm(q);
    if (q(1) < 0)
        q = q.*-1;
    end
    if (i < nSim)
        w = wFn(tOrb(i+1));
    end
end

% triad vector std test
siglist = siglist(:, any(siglist ~= 0, 1));

sig1 = sqrt(mean(siglist(1, :).^2));
sig2 = sqrt(mean(siglist(2, :).^2));
sigmag = mean(siglist(3, :));


% css fov lit test
nLitList = sum(yLitList == 1,1);
sensorCounts = [3;2;1];
lit = zeros(3,2);
for j = 1:3
    mask = nLitList == sensorCounts(j);
    lit(j,1) = nnz(mask);
    lit(j,2) = sum(thetaErrList(mask));
end


end