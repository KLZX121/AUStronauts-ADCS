%% plot_saved_mc.m
% Recreates the Monte Carlo plots from a previously saved mc_*.mat file.
% This script does NOT rerun the Monte Carlo simulations.

clear; close all; clc;

%% Select saved Monte Carlo results
[file, path] = uigetfile('../results/mc_*.mat', 'Select saved Monte Carlo results');
if isequal(file, 0)
    disp('No file selected.');
    return
end
load(fullfile(path, file));

hasB3 = exist('B3', 'var') && ~isempty(B3);
if ~hasB3
    B3 = zeros(size(A3));
end
%% Plot settings
set(0, 'DefaultLegendLocation', 'eastoutside')
set(0, 'DefaultLineLineWidth', 1.4)
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultTextFontSize', 12)

% Convergence close-up window (minutes). Edit these to suit.
CLOSEUP_TMIN = 0;
CLOSEUP_TMAX = 10;
BIAS_CLOSEUP_TMIN = 0;     % min
BIAS_CLOSEUP_TMAX = 20;    % min

% Recreate these from saved data if not present in the MAT file.
if ~exist('COMPARE_JOSEPH', 'var')
    COMPARE_JOSEPH = any(arrayfun(@(r) r.joseph, runs));
end
if ~exist('N_MC', 'var')
    N_MC = size(TH, 1);
end
if ~exist('labels', 'var') || isempty(labels)
    labels = strings(1, numel(runs));
    for k = 1:numel(runs)
        labels(k) = string(runs(k).filter) + " / " + string(runs(k).sunMode);
        if runs(k).joseph
            labels(k) = labels(k) + " / Joseph";
        end
    end
end

nRuns = size(TH, 3);	
tMin = sim.t ./ 60;
col = lines(2);                         % colour = filter
runCol = @(k) col((string(runs(k).filter) == "RIEKF") + 1, :);
runLS = @(k) LineStyleFor(runs(k), COMPARE_JOSEPH);

%% MC Attitude Error: full run
figure('Name', 'MC Attitude Error')
hold on
for k = 1:nRuns
    PlotBand(tMin, rad2deg(double(TH(:, :, k))), rad2deg(double(A3(:, :, k))), runCol(k), runLS(k), labels(k))
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('\theta error (deg)')
title(sprintf('Angular Error: RMS and 5-95%% band (%d trials, init: %s)', N_MC, cfg0.initAtt))
legend()
ShadeEclipse(tMin, sim.inEclipse)


%% MC Attitude Error: convergence close-up
figure('Name', 'MC Attitude Error - Convergence Close-up')
hold on
for k = 1:nRuns
    PlotBand(tMin, rad2deg(double(TH(:, :, k))), rad2deg(double(A3(:, :, k))), runCol(k), runLS(k), labels(k))
end
ShadeEclipse(tMin, sim.inEclipse)
hold off
yscale('log')
xlim([CLOSEUP_TMIN CLOSEUP_TMAX])
grid on
xlabel('t (min)')
ylabel('\theta error (deg)')
title(sprintf('Angular Error Convergence Close-up (%.1f-%.1f min)', CLOSEUP_TMIN, CLOSEUP_TMAX))
legend()

%% MC Gyro Bias Error
figure('Name', 'MC Gyro Bias Error')
hold on
for k = 1:nRuns
    PlotBand(tMin, rad2deg(double(BE(:, :, k))).*3600, rad2deg(double(B3(:, :, k))), runCol(k), runLS(k), labels(k))
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('|\Delta\beta| (deg/h)')
title('Gyro Bias Error Magnitude: RMS and 5-95% band')
legend()
ShadeEclipse(tMin, sim.inEclipse)

%% MC Gyro Bias Error (convergence close-up)

figure('Name', 'MC Gyro Bias Error - Convergence Close-up')
hold on

for k = 1:nRuns
    PlotBand(tMin, rad2deg(double(BE(:, :, k))).*3600, ...
        rad2deg(double(B3(:, :, k))), runCol(k), runLS(k), labels(k))
end

hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('|\Delta\beta| (deg/h)')
title('Gyro Bias Error: Convergence Close-up')

xlim([BIAS_CLOSEUP_TMIN BIAS_CLOSEUP_TMAX])
legend()
ShadeEclipse(tMin, sim.inEclipse)

%% NEES: median + band, total / attitude / bias
%{
neesSets = {NE, NA, NB};
neesNames = ["NEES (6 DOF)"; "Attitude NEES (3 DOF)"; "Bias NEES (3 DOF)"];
neesMedRef = [5.348, 2.366, 2.366];   % median chi-squared reference, 6/3/3 DOF
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
    yline(neesMedRef(p), 'k:', sprintf('chi^2 median (%.2f)', neesMedRef(p)), 'HandleVisibility', 'off')
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

%% NEES: mean over trials
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
    yline(neesMeanRef(p), 'k:', sprintf('Expected (%d)', neesMeanRef(p)), 'HandleVisibility', 'off')
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

%% Per-trial NEES traces
figure('Name', 'MC NEES (per trial)')
tl = tiledlayout('flow');
tl.Title.String = 'NEES (6 DOF): every trial (grey) and median (colour)';
for k = 1:nRuns
    nexttile
    hold on
    ShadeEclipse(tMin, sim.inEclipse)
    plot(tMin, double(NE(:, :, k))', 'Color', [0 0 0 0.15], 'LineWidth', 0.5, 'HandleVisibility', 'off')
    plot(tMin, median(double(NE(:, :, k)), 1, 'omitnan'), 'Color', runCol(k), 'LineWidth', 2)
    yline(6, 'k:', 'HandleVisibility', 'off')
    hold off
    yscale('log')
    grid on
    title(labels(k))
    xlabel('t (min)')
    ylabel('NEES')
end

%% Covariance health

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

%% Per-trial distributions
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
CONV_DEG = 5;
ylabel(sprintf('time to < %.0f deg (min)', CONV_DEG))
title(sprintf('Convergence time (%.0f%%-%.0f%% converged)', min(summary.ConvRatePct), max(summary.ConvRatePct)))
%}
%% Initial angular errors
initialErrors = nan(N_MC, nRuns);

for k = 1:nRuns
    for m = 1:N_MC
        if ~isnan(initIdx(m,k))
            initialErrors(m,k) = ...
                rad2deg(double(TH(m,initIdx(m,k),k)));
        end
    end
end

fprintf('Initial attitude error (deg):\n')
disp(array2table(initialErrors, 'VariableNames', ...
    matlab.lang.makeValidName(cellstr(labels))))

fprintf('RMS initial error (deg):\n')
disp(rms(initialErrors, 1, 'omitnan'))

fprintf('95th percentile initial error (deg):\n')
disp(prctile(initialErrors, 95, 1))

%% Local functions
function PlotBand(t, X, att3, c, ls, name)
% Plot rms with 5th-95th percentile band across Monte Carlo trials.
X = double(X);
med = rms(X, 1, 'omitnan');
lo = prctile(X, 5, 1);
hi = prctile(X, 95, 1);
ok = isfinite(med) & isfinite(lo) & isfinite(hi) & lo > 0;
fill([t(ok) fliplr(t(ok))], [lo(ok) fliplr(hi(ok))], c, ...
    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(t, med, ls, 'Color', c, 'DisplayName', name)

% plot sig-3
%plot(t, rms(att3, 1, 'omitnan'), ls, 'Color', 'k', 'DisplayName', name + ' 3\sigma')
end


function ls = LineStyleFor(run, compareJoseph)
% Line style by sun mode (and Joseph form when comparing).
if string(run.sunMode) == "vector"
    ls = '-';
else
    ls = ':';
end
if compareJoseph && run.joseph
    if string(run.sunMode) == "vector"
        ls = ':';
    else
        ls = '-.';
    end
end
end


function ShadeEclipse(tMin, inEclipse)
% Shade eclipse intervals on current axes.
dE = diff([0 double(inEclipse(:))' 0]);
s = find(dE == 1);
e = find(dE == -1) - 1;
for j = 1:numel(s)
    if s(j) <= numel(tMin) && e(j) <= numel(tMin)
        xregion(tMin(s(j)), tMin(e(j)), 'FaceColor', [0.5 0.5 0.5], ...
            'FaceAlpha', 0.15, 'HandleVisibility', 'off');
    end
end
end

