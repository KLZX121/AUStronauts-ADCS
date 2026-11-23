function out = RunSim(sim, cfg)
%RunSim Runs one 6-state attitude filter over pre-computed simulation data
%   out = RunSim(sim, cfg)
%
%   Inputs
%   sim             (struct)    Pre-computed truth + environment (run4sims.m)
%   .t              (1, n)      Time (s)
%   .dt             (1, 1)      Timestep (s)
%   .q              (4, n)      True attitude quaternion
%   .w              (3, n)      True body angular velocity (rad/s)
%   .bRef           (3, n)      Reference magnetic field (ECI, T)
%   .uSRef          (3, n)      Reference unit sun vector (ECI)
%   .inEclipse      (1, n)      Eclipse flag
%
%   cfg             (struct)    Run configuration
%   .filter         string      "MEKF" (MEKFNominal) or "RIEKF"
%   .sunMode        string      "scalar" (raw CSS) or "vector" (sun vector)
%   .seed           (1, 1)      rng seed (same for all runs => same noise)
%   .tInit          (1, 1)      Earliest filter initialisation time (s)
%
%   Optional (defaults in SetDefaults below)
%   .initAtt        string      Initial attitude estimate:
%                               "triad"    TRIAD from mag + sun (good)
%                               "identity" q = [1 0 0 0] (bad, arbitrary
%                                          error up to 180 deg)
%                               "random"   truth rotated by initErrDeg about
%                                          a random axis
%   .initErrDeg     (1, 1)      Error angle for "random" (deg)
%   .sigmaAttInit   (1, 1)      1-sigma initial attitude uncertainty for
%                               "identity" / "random" (rad, isotropic)
%   .sinMin         (1, 1)      Min |sin| of mag-sun angle to initialise
%                               (TRIAD only)
%   .sigmaSun       (1, 1)      Sun vector 1-sigma (rad), vector mode
%   .minSunSensors  (1, 1)      Min lit CSS for valid sun vector
%   .useJoseph      (1, 1)      Joseph-form (symmetrised) covariance update
%   .gyroFromQuat   (1, 1)      Feed the gyro the body rate implied by
%                               consecutive truth quaternions (ImpliedRate)
%                               instead of the trapezoid average of w(i-1)
%                               and w(i). Removes the sim-side integration
%                               mismatch (coning-type error) so that the
%                               filter's propagation is exactly consistent
%                               with the truth.
%   .noiseFree      (1, 1)      Disable all simulated sensor/gyro noise and
%                               bias random walk. The filter still uses its
%                               nominal Q and R, so this checks the filter
%                               and simulation for systematic errors.
%
%   Outputs
%   out             (struct)    Results for plotting / comparison

cfg = SetDefaults(cfg);

rng(cfg.seed)

% random initial attitude error axis, drawn up front so the noise stream
% seen by the sensors is identical for every filter/sun-mode in a trial
axisRand = randn(3, 1);
axisRand = axisRand/norm(axisRand);
aRand = deg2rad(cfg.initErrDeg);
dqRand = [cos(aRand/2); sin(aRand/2)*axisRand];

n = size(sim.q, 2);
dt = sim.dt;
t = sim.t;

%% sensor models (identical for every run)

% magnetometer
calM.bias = zeros(3, 1);
calM.D = zeros(3, 3);
calM.O = eye(3);
calM.sigma = (120/3)*1e-9;
mag = MagModel(calM);

% coarse sun sensors
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

% gyro (ideal scale factor / misalignment, as both filters are 6-state)
calG.bias = deg2rad([0.1; 0.1; 0.1])./3600;
calG.s = zeros(3, 1);
calG.kU = zeros(3, 1);
calG.kL = zeros(3, 1);
calG.sigG = repmat(sqrt(10)*1e-7, 3, 1);
calG.sigB = repmat(sqrt(10)*1e-10, 3, 1);
gyro = GyroModel(calG);

% nominal process noise for the filter, taken before any noise is disabled
QNom = diag([calG.sigG.^2; calG.sigB.^2]);
if cfg.noiseFree
    gyro.sigG = zeros(3, 1);
    gyro.sigB = zeros(3, 1);
end
genNoise = ~cfg.noiseFree;

%% filter

switch cfg.filter
    case "MEKF"
        f = MEKFNominal;
    case "RIEKF"
        f = RIEKF;
    otherwise
        error("RunSim:badFilter", "Unknown filter '%s'", cfg.filter)
end
f.sunMode = cfg.sunMode;
f.sigmaSun = cfg.sigmaSun;
f.minSunSensors = cfg.minSunSensors;
f.useJoseph = cfg.useJoseph;
f.genNoise = genNoise;

%% results
out.label = sprintf("%s / %s", cfg.filter, cfg.sunMode);
out.cfg = cfg;
out.t = t;
out.qEst = zeros(4, n);
out.thetaErr = nan(1, n);
out.xErr = zeros(f.nStates, n);
out.xSigma3 = zeros(f.nStates, n);
out.nees = nan(1, n);           % 6-state NEES
out.neesAtt = nan(1, n);        % attitude block NEES (3 DOF, marginal P)
out.neesBias = nan(1, n);       % bias block NEES (3 DOF, marginal P)
out.minDiagP = nan(1, n);       % min diag(P); <= 0 means P lost definiteness
out.minEigCorr = nan(1, n);     % min eig of correlation matrix of P (-1 if diag <= 0)
out.asymP = nan(1, n);          % ||P - P'||_F / ||P||_F
out.nSunUsed = zeros(1, n);
out.initIdx = NaN;
out.plotFn = str2func(class(f) + ".Plot");

wPrev = sim.w(:, 1);

if cfg.gyroFromQuat
    if isfield(sim, 'wImp')
        wImp = sim.wImp;
    else
        wImp = ImpliedRate(sim.q, dt);
    end
end

for i = 1:n
    q = sim.q(:, i);
    w = sim.w(:, i);

    env.bRef = sim.bRef(:, i);
    env.uSRef = sim.uSRef(:, i);
    env.T = 60;
    env.inEclipse = sim.inEclipse(i);

    %%% sensors (body frame)
    bMag = mag.Measurement(q, env.bRef, genNoise);
    uBMag = bMag/norm(bMag);

    ySMeas = css.Measurement(q, env, genNoise);
    uSMeas = css.CalcSunVec(ySMeas, env.T);

    if cfg.gyroFromQuat
        wAvg = wImp(:, i);
    else
        wAvg = 0.5*(wPrev + w);
    end
    wGyro = gyro.Measurement(wAvg, dt);
    wPrev = w;

    if ~f.initialised
        switch cfg.initAtt
            case "triad"
                uBRef = env.bRef./norm(env.bRef);
                [ATRIAD, PTRIAD] = TRIAD([uBMag, uSMeas], [uBRef, env.uSRef], mag.sigma/norm(bMag), 0.1);
                qEst = DCMToQ(ATRIAD);

                sinBS = norm(cross(uBMag, uSMeas));
                doInit = sinBS > cfg.sinMin && t(i) >= cfg.tInit;

                if cfg.filter == "RIEKF"
                    % inertial-frame attitude covariance
                    PAtt = ATRIAD'*PTRIAD*ATRIAD;
                else
                    PAtt = PTRIAD;
                end

            case "identity"
                qEst = [1; 0; 0; 0];
                doInit = t(i) >= cfg.tInit;
                PAtt = cfg.sigmaAttInit^2*eye(3);

            case "random"
                qEst = QProd(dqRand, q);
                if cfg.filter == "RIEKF"
                    qEst = RIEKF.HamProd(q, dqRand);
                end
                doInit = t(i) >= cfg.tInit;
                PAtt = cfg.sigmaAttInit^2*eye(3);
        end
        if qEst(1) < 0
            qEst = -qEst;
        end

        if doInit
            out.initIdx = i;

            x0 = zeros(6, 1);
            P0 = blkdiag(PAtt, (0.2*pi/(3600*180))^2.*eye(3));
            f.Initialise(x0, qEst, P0, QNom, gyro, mag, css);
        end
    else
        f.Propagate(dt, wGyro);
        f.Update(q, env);

        qEst = f.qEst;
        if qEst(1) < 0
            qEst = -qEst;
        end
        out.nSunUsed(i) = f.nVecMeas - 1;   % 1 if sun vector used (vector mode)
    end

    [out.thetaErr(i), ~] = QAttErr(q, qEst);
    out.qEst(:, i) = qEst;

    if f.initialised
        xErr = f.CalcError(q, gyro.bias);
        out.xErr(:, i) = xErr;
        out.xSigma3(:, i) = f.CalcSigma3();
        % NEES (6-state). Only meaningful once the attitude error is small:
        % the MEKF error vector is 2*tan(theta/2) and blows up near 180 deg
        out.nees(i) = xErr'*(f.P\xErr);

        % block NEES from the marginal covariances
        iA = 1:3;
        iB = 4:6;
        out.neesAtt(i) = xErr(iA)'*(f.P(iA, iA)\xErr(iA));
        out.neesBias(i) = xErr(iB)'*(f.P(iB, iB)\xErr(iB));

        % covariance health
        if cfg.logCov
            Pf = f.P;
            dP = diag(Pf);
            out.minDiagP(i) = min(dP);
            if all(dP > 0)
                sd = sqrt(dP);
                out.minEigCorr(i) = min(eig((Pf./(sd*sd') + (Pf./(sd*sd'))')/2));
            else
                out.minEigCorr(i) = -1;
            end
            out.asymP(i) = norm(Pf - Pf', 'fro')/norm(Pf, 'fro');
        end
    end

    gyro.PropagateBias(dt);
end

if cfg.sunMode == "scalar"
    out.nSunUsed = [];
end
end

function cfg = SetDefaults(cfg)
%SetDefaults Fills in optional cfg fields
def.initAtt = "triad";
def.initErrDeg = 90;
def.sigmaAttInit = deg2rad(60);
def.sinMin = 0.5;
def.sigmaSun = 0.1;
def.minSunSensors = 3;
def.tInit = 0;
def.useJoseph = false;
def.noiseFree = false;
def.gyroFromQuat = false;
def.logCov = true;

fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k})
        cfg.(fn{k}) = def.(fn{k});
    end
end
cfg.initAtt = string(cfg.initAtt);
end