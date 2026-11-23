%% run4simsMC: Monte Carlo MEKF vs RIEKF, sun vector vs raw coarse sun sensors
%
%   Run  Filter         Sun measurement
%   1    MEKFNominal    sun vector (CalcSunVec)
%   2    RIEKF          sun vector (CalcSunVec)
%   3    MEKFNominal    raw CSS outputs
%   4    RIEKF          raw CSS outputs
%
% With COMPARE_JOSEPH = true every configuration is also run with the
% Joseph-form (symmetrised) covariance update, giving 8 runs per trial.
%
% Each Monte Carlo trial m uses seed SEED0+m-1 for all runs, so within a
% trial the runs see the same sensor noise (paired comparison). The truth
% trajectory is fixed; only noise (and the random initial error axis when
% initAtt = "random") varies between trials.
%
% Initial attitude options (cfg0.initAtt):
%   "triad"     TRIAD estimate from mag + sun vector (good initial guess)
%   "identity"  q = [1 0 0 0], arbitrary error up to 180 deg
%   "random"    truth rotated by cfg0.initErrDeg about a random axis
%
% Diagnostics
%   - NEES: total (6 DOF), attitude block (3 DOF) and bias block (3 DOF),
%     as median + 5-95 % band, mean, and per-trial traces
%   - Covariance health: fraction of trials with non-positive-definite P,
%     median min-eigenvalue of the correlation matrix of P, max asymmetry
%   - cfg0.noiseFree = true: all simulated noise off (filter keeps its
%     nominal Q and R). Errors should go to ~0 with no growth; any
%     systematic error left is a simulation/filter mismatch, not tuning.
%     NEES is not meaningful in this mode.

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
if isempty(which('MagModel'))
    addpath(genpath('../sensors'));
end
if isempty(which('IGRFECI'))
    addpath(genpath('../environment'));
end

%% settings
SIM_LENGTH = 0;         % samples to use (0 = all)
N_MC = 1;              % number of Monte Carlo trials
SEED0 = 1;              % trial m uses seed SEED0 + m - 1
PAR_WORKERS = 0;        % 0 = serial, Inf = all workers (Parallel Computing Toolbox)
CONV_DEG = 5;           % "converged" = error stays below this after initialisation (deg)
COMPARE_JOSEPH = false; % true: run every config with short form AND Joseph form
SAVE_RESULTS = true;

cfg0.tInit = 0;                 % s, earliest initialisation time
cfg0.initAtt = "random";         % "triad" | "identity" | "random"
cfg0.initErrDeg = 60;           % deg, for "random"
cfg0.sigmaAttInit = deg2rad(cfg0.initErrDeg);% rad, initial attitude 1-sigma for "identity"/"random"
cfg0.sinMin = 0.5;              % min sin(mag-sun angle) to initialise (TRIAD only)
cfg0.sigmaSun = 0.1;            % rad, sun vector 1-sigma (vector mode)
cfg0.minSunSensors = 3;         % min lit CSS for a valid sun vector
cfg0.useJoseph = false;         % Joseph covariance update (ignored if COMPARE_JOSEPH)
cfg0.noiseFree = false;         % disable all simulated noise (see header)
cfg0.gyroFromQuat = true;      % gyro input = rate implied by truth quaternions (exact)
cfg0.logCov = false;             % log covariance health (small per-step cost)

%% read data
xData = jsondecode(readlines('../results/xdata.json', 'EmptyLineRule','skip'))';
jDData = jsondecode(readlines('../results/jd.json', 'EmptyLineRule', 'skip'));
eclipseData = jsondecode(readlines('../results/eclipse.json', 'EmptyLineRule', 'skip'));

if SIM_LENGTH > 0
    xData = xData(:, 1:SIM_LENGTH);
    jDData = jDData(1:SIM_LENGTH);
    eclipseData = eclipseData(1:SIM_LENGTH);
end

nSim = length(jDData);
tf = (jDData(end)-jDData(1))*86400;
dt = tf/(nSim-1);

d.iR = 1:3;
d.iQ = 7:10;
d.iWSat = 11:13;

%% pre-compute truth and environment once (shared by all runs and trials)
sim.t = 0:dt:tf;
sim.dt = dt;
sim.q = xData(d.iQ, :);
sim.w = xData(d.iWSat, :);
sim.bRef = zeros(3, nSim);
sim.uSRef = zeros(3, nSim);
sim.inEclipse = false(1, nSim);

for i = 1:nSim
    r = xData(d.iR, i);
    jD = jDData(i);

    sim.bRef(:, i) = IGRFECI(r, jD);
    [uS, rs] = SunV1(jD, r);
    sim.uSRef(:, i) = uS;
    sim.inEclipse(i) = CheckEclipse(r, uS.*rs);
end

% exact body rate between consecutive truth quaternions (shared by all runs)
if cfg0.gyroFromQuat
    sim.wImp = ImpliedRate(sim.q, dt);
end

%% configurations
base = struct('filter', {"MEKF", "RIEKF", "MEKF", "RIEKF"}, ...
              'sunMode', {"vector", "vector", "scalar", "scalar"});
if COMPARE_JOSEPH
    josephList = [false true];
else
    josephList = cfg0.useJoseph;
end

runs = struct('filter', {}, 'sunMode', {}, 'joseph', {});
for j = josephList
    for b = 1:numel(base)
        runs(end+1) = struct('filter', base(b).filter, ...
            'sunMode', base(b).sunMode, 'joseph', j); %#ok<SAGROW>
    end
end
nRuns = numel(runs);

labels = strings(1, nRuns);
for k = 1:nRuns
    labels(k) = runs(k).filter + " / " + runs(k).sunMode;
    if runs(k).joseph
        labels(k) = labels(k) + " / Joseph";
    end
end

%% Monte Carlo
fprintf("Monte Carlo: %d trials x %d runs, initAtt = %s, noiseFree = %d, gyroFromQuat = %d\n", ...
    N_MC, nRuns, cfg0.initAtt, cfg0.noiseFree, cfg0.gyroFromQuat)
tMC = tic;

mc = cell(N_MC, 1);
parfor (m = 1:N_MC, PAR_WORKERS)
    trial = cell(1, nRuns);
    for k = 1:nRuns
        cfg = cfg0;
        cfg.filter = runs(k).filter;
        cfg.sunMode = runs(k).sunMode;
        cfg.useJoseph = runs(k).joseph;
        cfg.seed = SEED0 + m - 1;

        r = RunSim(sim, cfg);

        % keep only what the statistics need (single to limit memory)
        trial{k} = struct( ...
            'thetaErr', single(r.thetaErr), ...
            'bErr', single(vecnorm(r.xErr(4:6, :))), ...
            'att3', single(vecnorm(r.xSigma3(1:3, :))), ...
            'bias3', single(vecnorm(r.xSigma3(4:6, :))), ...
            'nees', single(r.nees), ...
            'neesAtt', single(r.neesAtt), ...
            'neesBias', single(r.neesBias), ...
            'minDiagP', single(r.minDiagP), ...
            'minEigCorr', single(r.minEigCorr), ...
            'asymP', single(r.asymP), ...
            'initIdx', r.initIdx);
    end
    mc{m} = trial;
    fprintf("  trial %d/%d done\n", m, N_MC)
end
fprintf("Monte Carlo completed in %.1f s\n", toc(tMC))

%% stack results: (trial x time x run)
sz = [N_MC, nSim, nRuns];
TH = nan(sz, 'single');     % angular error (rad)
BE = nan(sz, 'single');     % |bias error| (rad/s)
A3 = nan(sz, 'single');     % ||attitude 3-sigma|| (rad)
B3 = nan(sz, 'single');
NE = nan(sz, 'single');     % NEES, 6 DOF
NA = nan(sz, 'single');     % NEES, attitude block, 3 DOF
NB = nan(sz, 'single');     % NEES, bias block, 3 DOF
MD = nan(sz, 'single');     % min diag(P)
ME = nan(sz, 'single');     % min eig of correlation matrix of P
AS = nan(sz, 'single');     % asymmetry of P
initIdx = nan(N_MC, nRuns);
for m = 1:N_MC
    for k = 1:nRuns
        s = mc{m}{k};
        TH(m, :, k) = s.thetaErr;
        BE(m, :, k) = s.bErr;
        A3(m, :, k) = s.att3;
        B3(m, :, k) = s.bias3;
        NE(m, :, k) = s.nees;
        NA(m, :, k) = s.neesAtt;
        NB(m, :, k) = s.neesBias;
        MD(m, :, k) = s.minDiagP;
        ME(m, :, k) = s.minEigCorr;
        AS(m, :, k) = s.asymP;
        initIdx(m, k) = s.initIdx;
    end
end
clear mc s

%% per-trial metrics
convTime = nan(N_MC, nRuns);    % s from initialisation until error stays < CONV_DEG
ssRMS = nan(N_MC, nRuns);       % deg, RMS over second half after initialisation
ssNEES = nan(N_MC, nRuns);      % mean NEES (6 DOF) over second half after init
ssNEESAtt = nan(N_MC, nRuns);   % same, attitude block
ssNEESBias = nan(N_MC, nRuns);  % same, bias block
lateNEES = nan(N_MC, nRuns);    % mean NEES (6 DOF) over final 10 % of the run
notPD = false(N_MC, nRuns);     % P lost positive definiteness at any time
maxAsym = nan(N_MC, nRuns);     % worst asymmetry of P over the run

late = round(0.9*nSim):nSim;

for k = 1:nRuns
    for m = 1:N_MC
        i0 = initIdx(m, k);
        if isnan(i0)
            continue
        end
        e = rad2deg(double(TH(m, i0:end, k)));

        bad = find(~(e < CONV_DEG));
        if isempty(bad)
            convTime(m, k) = 0;
        elseif bad(end) < numel(e)
            convTime(m, k) = sim.t(i0 + bad(end)) - sim.t(i0);
        end

        half = round(numel(e)/2):numel(e);
        iHalf = i0 + half - 1;
        ssRMS(m, k) = sqrt(mean(e(half).^2, 'omitnan'));
        ssNEES(m, k) = mean(double(NE(m, iHalf, k)), 'omitnan');
        ssNEESAtt(m, k) = mean(double(NA(m, iHalf, k)), 'omitnan');
        ssNEESBias(m, k) = mean(double(NB(m, iHalf, k)), 'omitnan');
        lateNEES(m, k) = mean(double(NE(m, late, k)), 'omitnan');

        notPD(m, k) = any(ME(m, :, k) <= 0) || any(MD(m, :, k) <= 0);
        maxAsym(m, k) = max(AS(m, :, k), [], 'omitnan');
    end
end

%% summary table
vars = ["ConvRatePct" "ConvMedMin" "ConvP95Min" "SSRmsMedDeg" "SSRmsP95Deg" ...
        "MeanNEES" "MedNEES" "MedNEESAtt" "MedNEESBias" "LateMedNEES" ...
        "NotPDTrials" "MaxAsym" "InitMedMin" "NTrials"];
summary = table('Size', [nRuns numel(vars)], ...
    'VariableTypes', repmat("double", 1, numel(vars)), ...
    'VariableNames', vars, 'RowNames', labels);
for k = 1:nRuns
    ct = convTime(:, k);
    summary.ConvRatePct(k) = 100*mean(isfinite(ct));
    summary.ConvMedMin(k) = Pct(ct, 50)/60;
    summary.ConvP95Min(k) = Pct(ct, 95)/60;
    summary.SSRmsMedDeg(k) = Pct(ssRMS(:, k), 50);
    summary.SSRmsP95Deg(k) = Pct(ssRMS(:, k), 95);
    summary.MeanNEES(k) = mean(ssNEES(:, k), 'omitnan');
    summary.MedNEES(k) = Pct(ssNEES(:, k), 50);
    summary.MedNEESAtt(k) = Pct(ssNEESAtt(:, k), 50);
    summary.MedNEESBias(k) = Pct(ssNEESBias(:, k), 50);
    summary.LateMedNEES(k) = Pct(lateNEES(:, k), 50);
    summary.NotPDTrials(k) = sum(notPD(:, k));
    summary.MaxAsym(k) = max(maxAsym(:, k), [], 'omitnan');
    okI = isfinite(initIdx(:, k));
    summary.InitMedMin(k) = Pct(sim.t(initIdx(okI, k)), 50)/60;
    summary.NTrials(k) = N_MC;
end
fprintf("\nConvergence threshold: %.1f deg\n", CONV_DEG)
fprintf("Expected NEES for a consistent filter: 6 (total), 3 (each block)\n")
fprintf("NEES columns are per-trial second-half means, then mean/median over trials.\n")
fprintf("LateMedNEES: median over trials of the mean NEES in the final 10%% of the run.\n")
if cfg0.noiseFree
    fprintf("NOISE-FREE mode: NEES is not meaningful; look at the error plots.\n")
end
disp(summary)

if SAVE_RESULTS
    tag = "";
    if cfg0.noiseFree
        tag = tag + "_nf";
    end
    if cfg0.gyroFromQuat
        tag = tag + "_gq";
    end
    if COMPARE_JOSEPH
        tag = tag + "_joseph";
    elseif cfg0.useJoseph
        tag = tag + "_j";
    end
    save(sprintf('../results/mc_%s%s.mat', cfg0.initAtt, tag), ...
        'TH', 'BE', 'A3', 'B3', 'NE', 'NA', 'NB', 'MD', 'ME', 'AS', 'initIdx', ...
        'convTime', 'ssRMS', 'ssNEES', 'ssNEESAtt', 'ssNEESBias', 'lateNEES', ...
        'notPD', 'maxAsym', 'summary', 'cfg0', 'runs', 'labels', ...
        'N_MC', 'SEED0', 'COMPARE_JOSEPH', 'sim', '-v7.3')
end

%% plots
tMin = sim.t./60;
col = lines(2);                         % colour = filter
runCol = @(k) col((runs(k).filter == "RIEKF") + 1, :);
runLS = @(k) LineStyleFor(runs(k), COMPARE_JOSEPH);

% median + 5-95 percentile bands
figure('Name', 'MC Attitude Error')
hold on
for k = 1:nRuns
    PlotBand(tMin, rad2deg(double(TH(:, :, k))), runCol(k), runLS(k), labels(k))
    sig3Deg = rms(rad2deg(double(A3(:, :, k))), 1, 'omitnan');
    plot(tMin, sig3Deg, runLS(k), 'Color', 'k', 'DisplayName', labels(k) + " 3\sigma")
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('\theta error (deg)')
title(sprintf('Angular Error: median and 5-95%% band (%d trials, init: %s)', N_MC, cfg0.initAtt))
legend()
PlotEclipsePatch(eclipseData, tMin)

figure('Name', 'MC Gyro Bias Error')
hold on
for k = 1:nRuns
    PlotBand(tMin, rad2deg(double(BE(:, :, k))).*3600, runCol(k), runLS(k), labels(k))
    b3DegHr = rms(rad2deg(double(B3(:, :, k))).*3600, 1, 'omitnan');
    plot(tMin, b3DegHr, runLS(k), 'Color', 'k', 'DisplayName', labels(k) + " 3\sigma")
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('|\Delta\beta| (deg/h)')
title('Gyro Bias Error Magnitude: median and 5-95% band')
legend()
PlotEclipsePatch(eclipseData, tMin)

% NEES: median + band, total / attitude / bias
neesSets = {NE, NA, NB};
neesNames = ["NEES (6 DOF)"; "Attitude NEES (3 DOF)"; "Bias NEES (3 DOF)"];
neesMedRef = [5.348, 2.366, 2.366];   % median of chi-squared with 6, 3, 3 DOF
neesMeanRef = [6, 3, 3];

figure('Name', 'MC NEES (median)')
tl = tiledlayout(3, 1);
tl.Title.String = 'NEES: median and 5-95% band over trials';
for p = 1:3
    nexttile
    hold on
    ShadeEclipse(tMin, sim.inEclipse)
    for k = 1:nRuns
        PlotBand(tMin, double(neesSets{p}(:, :, k)), runCol(k), runLS(k), labels(k))
    end
    yline(neesMedRef(p), 'k:', sprintf('chi^2 median (%.2f)', neesMedRef(p)), ...
        'HandleVisibility', 'off')
    hold off
    yscale('log')
    grid on
    ylabel(neesNames(p))
    if p == 1
        legend()
    end
    if p < 3
        xticklabels({})
    else
        xlabel('t (min)')
    end
end

figure('Name', 'MC NEES (mean)')
tl = tiledlayout(3, 1);
tl.Title.String = 'NEES: mean over trials (sensitive to outlier trials)';
for p = 1:3
    nexttile
    hold on
    ShadeEclipse(tMin, sim.inEclipse)
    for k = 1:nRuns
        plot(tMin, mean(double(neesSets{p}(:, :, k)), 1, 'omitnan'), runLS(k), ...
            'Color', runCol(k), 'DisplayName', labels(k))
    end
    yline(neesMeanRef(p), 'k:', sprintf('Expected (%d)', neesMeanRef(p)), ...
        'HandleVisibility', 'off')
    hold off
    yscale('log')
    grid on
    ylabel(neesNames(p))
    if p == 1
        legend()
    end
    if p < 3
        xticklabels({})
    else
        xlabel('t (min)')
    end
end

% per-trial NEES traces (is the late rise a few trials or all of them?)
figure('Name', 'MC NEES (per trial)')
tl = tiledlayout('flow');
tl.Title.String = 'NEES (6 DOF): every trial (grey) and median (colour)';
for k = 1:nRuns
    nexttile
    hold on
    ShadeEclipse(tMin, sim.inEclipse)
    plot(tMin, double(NE(:, :, k))', 'Color', [0 0 0 0.15], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off')
    plot(tMin, median(double(NE(:, :, k)), 1, 'omitnan'), 'Color', runCol(k), 'LineWidth', 2)
    yline(6, 'k:', 'HandleVisibility', 'off')
    hold off
    yscale('log')
    grid on
    title(labels(k))
    xlabel('t (min)')
    ylabel('NEES')
end

% covariance health
figure('Name', 'MC Covariance Health')
tl = tiledlayout(3, 1);
tl.Title.String = 'Covariance health';

nexttile
hold on
ShadeEclipse(tMin, sim.inEclipse)
for k = 1:nRuns
    nonPD = (double(ME(:, :, k)) <= 0) | (double(MD(:, :, k)) <= 0);
    plot(tMin, mean(nonPD, 1), runLS(k), 'Color', runCol(k), 'DisplayName', labels(k))
end
hold off
grid on
ylabel('frac. trials P not PD')
ylim([-0.05 1.05])
legend()
xticklabels({})

nexttile
hold on
ShadeEclipse(tMin, sim.inEclipse)
for k = 1:nRuns
    plot(tMin, median(double(ME(:, :, k)), 1, 'omitnan'), runLS(k), ...
        'Color', runCol(k), 'DisplayName', labels(k))
end
hold off
yscale('log')
grid on
ylabel('median min eig(corr P)')
xticklabels({})

nexttile
hold on
ShadeEclipse(tMin, sim.inEclipse)
for k = 1:nRuns
    plot(tMin, max(max(double(AS(:, :, k)), [], 1, 'omitnan'), 1e-20), runLS(k), ...
        'Color', runCol(k), 'DisplayName', labels(k))
end
hold off
yscale('log')
grid on
ylabel('max asym(P)')
xlabel('t (min)')

% per-trial distributions
figure('Name', 'MC Distributions')
tiledlayout(1, 2);

nexttile
g = repmat(labels, N_MC, 1);
boxchart(categorical(g(:), labels, labels), ssRMS(:))
yscale('log')
grid on
ylabel('steady-state RMS error (deg)')
title('Steady-state error (second half)')

nexttile
boxchart(categorical(g(:), labels, labels), convTime(:)./60)
grid on
ylabel(sprintf('time to < %.0f deg (min)', CONV_DEG))
title(sprintf('Convergence time (%.0f%%-%.0f%% converged)', ...
    min(summary.ConvRatePct), max(summary.ConvRatePct)))

%% local functions
function p = Pct(x, pr)
%Pct Percentile ignoring NaN; NaN if no valid samples
x = double(x(isfinite(x)));
if isempty(x)
    p = NaN;
else
    p = prctile(x, pr);
end
end

function PlotBand(t, X, c, ls, name)
%PlotBand Plots median of X (trials x time) with a 5-95 percentile band
X = double(X);
med = median(X, 1, 'omitnan');
lo = prctile(X, 5, 1);
hi = prctile(X, 95, 1);
ok = isfinite(med) & isfinite(lo) & isfinite(hi) & lo > 0;

fill([t(ok) fliplr(t(ok))], [lo(ok) fliplr(hi(ok))], c, ...
    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(t, med, ls, 'Color', c, 'DisplayName', name)
end

function ls = LineStyleFor(run, compareJoseph)
%LineStyleFor Line style by sun mode (and Joseph form when comparing)
if run.sunMode == "vector"
    ls = '-';
else
    ls = '--';
end
if compareJoseph && run.joseph
    if run.sunMode == "vector"
        ls = ':';
    else
        ls = '-.';
    end
end
end

function ShadeEclipse(tMin, inEclipse)
%ShadeEclipse Shades eclipse intervals on the current axes
dE = diff([0 double(inEclipse(:))' 0]);
s = find(dE == 1);
e = find(dE == -1) - 1;
for j = 1:numel(s)
    xregion(tMin(s(j)), tMin(e(j)), 'FaceColor', [0.5 0.5 0.5], ...
        'FaceAlpha', 0.15, 'HandleVisibility', 'off');
end
end