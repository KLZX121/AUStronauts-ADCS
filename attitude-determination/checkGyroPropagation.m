%% checkGyroPropagation: is the truth data consistent with the filter's propagation?
%
% Noise-free test of the simulation <-> filter interface. No sensors, no
% noise, no bias, no measurement updates. It answers:
%
%   1. Is the time base uniform? (the filter assumes dt = tf/(n-1))
%   2. If the true body rate w(i) is integrated with the filter's own
%      MEKF.PropQDisc, how far does the attitude drift from the truth
%      quaternion, for different ways of choosing the rate per step?
%        - trapezoid:       0.5*(w(i)   + w(i+1))   <- what main/RunSim feed
%        - previous sample: w(i)
%        - current sample:  w(i+1)
%   3. What constant-rate "gyro error" does that correspond to? The rate
%      implied by consecutive truth quaternions is compared with each
%      variant. Its mean (deg/h) is the bias-like error the filter would
%      see; compare it with the true gyro bias of 0.1 deg/h per axis.
%
% Reading the results
%   - Mean mismatch much smaller than 0.1 deg/h for some variant: that
%     variant is consistent with the truth and the interface is fine.
%   - Mean mismatch comparable to or larger than 0.1 deg/h for the variant
%     used in the sims: a systematic error is being absorbed into the bias
%     estimate (the late bias NEES rise).
%   - Large relative error between implied and applied rate: the
%     quaternion convention assumed here (q_new = q_old (x) dq, Hamilton)
%     does not match xdata.json, so the other numbers are meaningless.

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
SIM_LENGTH = 0;             % samples to use (0 = all)
BIAS_REF_DEGH = 0.1;        % deg/h, true gyro bias per axis used in the sims
SMOOTH_S = 60;              % s, moving-average window for the rate plots

%% read data
xData = jsondecode(readlines('../results/xdata.json', 'EmptyLineRule','skip'))';
jDData = jsondecode(readlines('../results/jd.json', 'EmptyLineRule', 'skip'));

if SIM_LENGTH > 0
    xData = xData(:, 1:SIM_LENGTH);
    jDData = jDData(1:SIM_LENGTH);
end

d.iR = 1:3;
d.iQ = 7:10;
d.iWSat = 11:13;

n = length(jDData);
tf = (jDData(end) - jDData(1))*86400;
dtNom = tf/(n-1);               % what the sims/filter use
t = (0:n-1)*dtNom;

q = xData(d.iQ, :);
w = xData(d.iWSat, :);

% eclipse flags for shading (same calls as run4simsMC.m)
inEclipse = false(1, n);
for i = 1:n
    r = xData(d.iR, i);
    [uS, rs] = SunV1(jDData(i), r);
    inEclipse(i) = CheckEclipse(r, uS.*rs);
end

%% 1. time base
dtVec = diff(jDData(:))'*86400;             % s, 1 x (n-1)
dtErr = dtVec - dtNom;
tOffset = [0 cumsum(dtErr)];                % s, accumulated offset vs uniform grid
jdRes = eps(jDData(1))*86400;               % s, resolution of a Julian date double

fprintf("\n=== Time base ===\n")
fprintf("Nominal dt (filter)          : %.9f s\n", dtNom)
fprintf("Per-step dt: min %.6f, max %.6f, std %.2e s\n", min(dtVec), max(dtVec), std(dtVec))
fprintf("Julian date resolution       : %.2e s (per-step dt jitter below this is quantisation)\n", jdRes)
fprintf("Max accumulated offset       : %.3e s\n", max(abs(tOffset)))
fprintf("Relative dt error (rms)      : %.2e\n", rms(dtErr)/dtNom)

%% 2. open-loop propagation of the true rate with the filter's own code
% rate used for the step i -> i+1; the sims compute wAvg at loop index i+1
% as 0.5*(w(i) + w(i+1))
variants = ["trapezoid (sim)"; "previous sample"; "current sample"];
wVar = {
    0.5*(w(:, 1:end-1) + w(:, 2:end));
    w(:, 1:end-1);
    w(:, 2:end);
};
nVar = numel(wVar);

thetaErr = zeros(nVar + 1, n);              % rad
for v = 1:nVar
    qEst = OpenLoop(q(:, 1), wVar{v}, dtNom);
    thetaErr(v, :) = AngErr(q, qEst);
end
% trapezoid with the true per-step dt instead of the nominal one
qEst = OpenLoop(q(:, 1), wVar{1}, dtVec);
thetaErr(nVar + 1, :) = AngErr(q, qEst);
varNames = [variants; "trapezoid, true dt"];

fprintf("\n=== Open-loop attitude drift (no measurements) ===\n")
tH = t(end)/3600;
for v = 1:nVar + 1
    fprintf("%-20s final %9.5f deg, max %9.5f deg, final/duration %8.4f deg/h\n", ...
        varNames(v), rad2deg(thetaErr(v, end)), rad2deg(max(thetaErr(v, :))), ...
        rad2deg(thetaErr(v, end))/tH)
end

%% 3. rate implied by consecutive truth quaternions vs applied rate
qConj = @(p) [p(1); -p(2:4)];
wImp = zeros(3, n-1);           % rad/s, body frame, using the filter's dt
for i = 1:n-1
    dq = RIEKF.HamProd(qConj(q(:, i)), q(:, i+1));
    dq = dq/norm(dq);
    wImp(:, i) = RIEKF.QLog(dq)/dtNom;
end

wTrap = wVar{1};
relErr = vecnorm(wImp - wTrap)./max(vecnorm(wTrap), eps);
fprintf("\n=== Implied rate (from q) vs trapezoid rate ===\n")
fprintf("Median relative difference: %.2e (should be << 1)\n", median(relErr))
if median(relErr) > 0.05
    warning("checkGyroPropagation:convention", ...
        ["Implied and applied rates disagree by >5%%: the assumed quaternion " ...
         "convention probably does not match xdata.json. Results below are not meaningful."])
end

segEcl = inEclipse(1:end-1);
tabVars = ["MeanX" "MeanY" "MeanZ" "MeanNorm" "RMS" "MeanNormSun" "MeanNormEcl" "ScalePPM"];
mismatch = table('Size', [nVar numel(tabVars)], ...
    'VariableTypes', repmat("double", 1, numel(tabVars)), ...
    'VariableNames', tabVars, 'RowNames', variants);
dwAll = cell(nVar, 1);
for v = 1:nVar
    dw = rad2deg(wImp - wVar{v})*3600;      % deg/h
    dwAll{v} = dw;
    m = mean(dw, 2);
    mismatch.MeanX(v) = m(1);
    mismatch.MeanY(v) = m(2);
    mismatch.MeanZ(v) = m(3);
    mismatch.MeanNorm(v) = norm(m);
    mismatch.RMS(v) = sqrt(mean(sum(dw.^2, 1)));
    mismatch.MeanNormSun(v) = norm(mean(dw(:, ~segEcl), 2));
    mismatch.MeanNormEcl(v) = norm(mean(dw(:, segEcl), 2));
    % least-squares scale between implied and applied rate
    mismatch.ScalePPM(v) = 1e6*(sum(wImp.*wVar{v}, 'all')/sum(wVar{v}.^2, 'all') - 1);
end
fprintf("Mismatch = implied rate - applied rate, in deg/h (true bias is %.2f deg/h per axis, %.3f deg/h norm)\n", ...
    BIAS_REF_DEGH, sqrt(3)*BIAS_REF_DEGH)
disp(mismatch)

%% plots
tMin = t./60;

figure('Name', 'Time Base')
tiledlayout(2, 1);
nexttile
plot(tMin(2:end), dtErr, '.', 'MarkerSize', 4)
hold on
yline([-jdRes jdRes], 'r--', 'HandleVisibility', 'off')
hold off
grid on
ylabel('dt - dt_{nom} (s)')
title('Per-step time error (red: Julian date resolution)')
xticklabels({})
nexttile
plot(tMin, tOffset)
grid on
xlabel('t (min)')
ylabel('accumulated offset (s)')

figure('Name', 'Open-loop Attitude Drift')
hold on
ShadeEclipse(tMin, inEclipse)
for v = 1:nVar + 1
    plot(tMin, rad2deg(thetaErr(v, :)), 'DisplayName', varNames(v))
end
hold off
yscale('log')
grid on
xlabel('t (min)')
ylabel('\theta error (deg)')
title('Open-loop propagation of true rate vs truth quaternion')
legend()

figure('Name', 'Rate Mismatch')
tiledlayout(4, 1);
axLbl = ["x"; "y"; "z"];
win = max(1, round(SMOOTH_S/dtNom));
for a = 1:3
    nexttile
    hold on
    ShadeEclipse(tMin, inEclipse)
    for v = 1:nVar
        plot(tMin(1:end-1), movmean(dwAll{v}(a, :), win), 'DisplayName', variants(v))
    end
    yline([-BIAS_REF_DEGH BIAS_REF_DEGH], 'k:', 'HandleVisibility', 'off')
    hold off
    grid on
    ylabel(sprintf('\\Delta\\omega_%s (deg/h)', axLbl(a)))
    if a == 1
        title(sprintf('Implied - applied rate, %d s moving mean (dotted: true bias)', SMOOTH_S))
        legend()
    end
    xticklabels({})
end
nexttile
hold on
ShadeEclipse(tMin, inEclipse)
plot(tMin, rad2deg(vecnorm(w)))
hold off
grid on
ylabel('|\omega| (deg/s)')
xlabel('t (min)')

%% local functions
function qEst = OpenLoop(q0, wMat, dt)
%OpenLoop Propagates q0 through the rates in wMat with MEKF.PropQDisc
%   dt may be a scalar or a vector with one value per step
nn = size(wMat, 2) + 1;
qEst = zeros(4, nn);
qEst(:, 1) = q0;
for i = 1:nn-1
    dti = dt(min(i, numel(dt)));
    wi = wMat(:, i);
    if norm(wi) < 1e-12
        qEst(:, i+1) = qEst(:, i);      % PropQDisc divides by |w|
    else
        qEst(:, i+1) = MEKF.PropQDisc(qEst(:, i), wi, dti);
    end
end
end

function th = AngErr(q, qEst)
%AngErr Angular error (rad) between two quaternion histories
n = size(q, 2);
th = zeros(1, n);
for i = 1:n
    th(i) = QAttErr(q(:, i), qEst(:, i));
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