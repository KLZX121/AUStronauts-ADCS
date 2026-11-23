%% run4sims: MEKF vs RIEKF, constructed sun vector vs raw coarse sun sensors
%
%   Run  Filter         Sun measurement          Invariant?
%   1    MEKFNominal    sun vector (CalcSunVec)  n/a
%   2    RIEKF          sun vector (CalcSunVec)  yes (vector measurement)
%   3    MEKFNominal    raw CSS outputs          n/a
%   4    RIEKF          raw CSS outputs          no  (scalar measurement)
%
% The magnetometer is a vector measurement so it is invariant in all runs.
% Every run uses the same rng seed and consumes the same number of random
% draws per step, so all four see the same noise realisation.

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
PLOT_EACH_RUN = false;  % per-run overview + detailed state error plots

cfg0.seed = 1;
cfg0.tInit = 0;             % s, earliest initialisation time
cfg0.sinMin = 0.5;          % min sin(mag-sun angle) to initialise
cfg0.sigmaSun = 0.1;        % rad, sun vector 1-sigma (vector mode)
cfg0.minSunSensors = 3;     % min lit CSS for a valid sun vector

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

%% pre-compute truth and environment once (shared by all runs)
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

%% configurations
runs = struct('filter', {"MEKF", "RIEKF", "MEKF", "RIEKF"}, ...
              'sunMode', {"vector", "vector", "scalar", "scalar"});

results = cell(1, numel(runs));
for k = 1:numel(runs)
    cfg = cfg0;
    cfg.filter = runs(k).filter;
    cfg.sunMode = runs(k).sunMode;

    fprintf("Run %d: %s / %s ...\n", k, cfg.filter, cfg.sunMode)
    tic
    results{k} = RunSim(sim, cfg);
    fprintf("   done in %.1f s, initialised at %.2f min\n", toc, ...
        sim.t(results{k}.initIdx)/60)
end

%% summary table
summary = table('Size', [numel(runs) 5], ...
    'VariableTypes', ["string" "double" "double" "double" "double"], ...
    'VariableNames', ["Run" "InitMin" "MedianDeg" "RMSDeg" "FinalDeg"]);
for k = 1:numel(runs)
    r = results{k};
    idx = r.initIdx:nSim;
    idx = idx(round(numel(idx)/2):end);   % second half after init
    e = rad2deg(r.thetaErr(idx));
    e = e(isfinite(e));

    summary.Run(k) = r.label;
    summary.InitMin(k) = sim.t(r.initIdx)/60;
    summary.MedianDeg(k) = median(e);
    summary.RMSDeg(k) = sqrt(mean(e.^2));
    summary.FinalDeg(k) = rad2deg(r.thetaErr(end));
end
disp(summary)

%% comparison plots
tMin = sim.t./60;
col = lines(2);                 % colour = filter
ls = containers.Map({'vector', 'scalar'}, {'-', '--'}); % style = sun mode

figure('Name', 'Attitude Error Comparison')
hold on
for k = 1:numel(runs)
    r = results{k};
    c = col((r.cfg.filter == "RIEKF") + 1, :);
    plot(tMin, rad2deg(r.thetaErr), ls(char(r.cfg.sunMode)), 'Color', c, ...
        'DisplayName', r.label)
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('\theta error (deg)')
title('Angular Error: MEKF vs RIEKF')
legend()
PlotEclipsePatch(eclipseData, tMin)

figure('Name', 'Gyro Bias Error Comparison')
hold on
for k = 1:numel(runs)
    r = results{k};
    c = col((r.cfg.filter == "RIEKF") + 1, :);
    bErr = rad2deg(vecnorm(r.xErr(4:6, :))).*3600;
    plot(tMin, bErr, ls(char(r.cfg.sunMode)), 'Color', c, 'DisplayName', r.label)
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('|\Delta\beta| (deg/h)')
title('Gyro Bias Error Magnitude')
legend()
PlotEclipsePatch(eclipseData, tMin)

figure('Name', 'Attitude 3-sigma Comparison')
hold on
for k = 1:numel(runs)
    r = results{k};
    c = col((r.cfg.filter == "RIEKF") + 1, :);
    % norm of the per-axis 3-sigma bounds (frame independent magnitude)
    s3 = rad2deg(vecnorm(r.xSigma3(1:3, :)));
    s3(s3 == 0) = NaN;
    plot(tMin, s3, ls(char(r.cfg.sunMode)), 'Color', c, 'DisplayName', r.label)
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('||3\sigma|| (deg)')
title('Attitude Covariance (3-sigma norm)')
legend()

if PLOT_EACH_RUN
    for k = 1:numel(runs)
        r = results{k};
        r.plotFn(tMin, r.xErr, r.xSigma3, 'overview')
        r.plotFn(tMin, r.xErr, r.xSigma3, 'detailed')
    end
end