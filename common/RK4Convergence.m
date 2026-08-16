function out = RK4Convergence(x0, xDotFn, d, nTests, nSigFig, dTest)
%RK4Convergence Runs a step size study for RK4
%   out = RK4Convergence(x0, xDotFn, d, nTests, nSigFig, dTest)
%   Use to determine optimal step size for RK4
%   Use out.hList(out.hStop) to get step size at which required relative 
%   tolerance is met. out.hStop will be 0 if tolerance is not met for all 
%   state values.Increase nTests if tolerance is not yet met. 
%   Plots will show approxErr, relErr, convergence ratio and observed 
%   order (p). Convergence ratio should approach 16 and p should approach 4
%   Note that values will diverge if approxErr approaches 0.
%
%   Inputs:
%   x0      (nx, 1)     Initial state values
%   xDotFn  (:, 1)      Functions for state derivatives
%   d       (struct)    Input struct for functions
%   nTests              Number of convergence steps (each test halves step
%                       size)
%   nSigFig             Result will be accurate to at least nSigFig, 
%                       calculated using the Scarborough criterion
%                       (optional - default 3)
%   dTest   (struct)    Test parameters (optional)
%           .t0         Initial time (default 0)
%           .tf         Final time  (default 1)
%           .h0         Initial step (default 1)
%
%   Outputs:
%   out     (struct)    Ouput struct
%           .hStop                      Step index at which relative error 
%                                       is below relative tolerance
%                                       (nSigFig)
%           .xList      (nx, nTests)    List of states for each test
%           .hList      (1, nTests)     List of steps
%           .approxErr  (nx, nTests-1)  Approximate error (x_i - x_(i-1))
%           .relErr     (nx, nTests-1)  Relative error (approxErr/x_i)

if (nargin < 5)
    nSigFig = 3;
end
if (nargin < 6)
    dTest.t0 = 0;
    dTest.tf = 1;
    dTest.h0 = 1;
end

% good implementation of rk4 in chapra pg. 638
% verify with ode45

out.xList = zeros(length(x0), nTests);
out.hList = zeros(1, nTests);
out.hList(1) = dTest.h0;

% reduction formula
hNew = @(h) 0.5*h;

for testI = 1:nTests
    % get step
    h = out.hList(testI);

    % run RK4 with step
    N = (dTest.tf-dTest.t0)/h;

    fprintf("Test %d: step size = %f steps = %d\n", testI, h, N)
    
    x = x0;
    
    tic
    for i = 1:N
        x = PropState(xDotFn, x, d, h);
    end
    out.xList(:, testI) = x;
    simTime = toc;
    fprintf("        time = %f s\n", simTime)

    % generate next step if not last step already
    if (testI < nTests)
        out.hList(testI+1) = hNew(h);
    end
end

% approx error, relative error: converge -> 0
out.approxErr = zeros(length(x), nTests-1);
out.relErr = zeros(size(out.approxErr));

% stopping criterion
% note relTol is in % while relErr is not.
% src - chapra pg. 104
relTolFn = @(n) 0.5*10^(-n);
relTol = relTolFn(nSigFig);
fprintf("\nLeast amount of accurate sig figs is %d\n", nSigFig)
fprintf("Relative tolerance is %e\n\n", relTol);

% step at which relTol is met
out.hStop = 0;

for i = 1:nTests-1
    out.approxErr(:, i) = abs(out.xList(:, i+1) - out.xList(:, i));
    out.relErr(:, i) = abs(out.approxErr(:, i) ./ out.xList(:, i+1));

    % find the maximum relative error in the state vector
    % ignore values that are less than machine epsilon as their errors are
    % unreliable
    findMaxRel = true;
    scanErrList = out.relErr(:, i);
    ignoredCount = 0;
    while (findMaxRel)
        [relMax, relI] = max(scanErrList);
        relMaxValue = out.xList(relI, i+1);

        if (abs(relMaxValue) > eps) 
            findMaxRel = false;
        else
            ignoredCount = ignoredCount + 1;
            scanErrList(relI) = 0;
        end
    end
    fprintf("Test %d: max rel err: %e (state %d) val = %e\n", i+1, relMax, relI, relMaxValue);
    fprintf("    Number of ignored values due to machine precision: %d\n", ignoredCount)

    % if the found relative error is below tolerance then set recommended
    % step size
    if ((out.hStop == 0) && (relMax < relTol))  
        out.hStop = i+1;
    end
end
if (~out.hStop)
    disp("No appropriate step was found - manually check values and plots")
else
    disp(' ')
    disp('Tolerance met at step ' + string(out.hStop))
    disp('Recommended step size is ' + string(out.hList(out.hStop)))
    disp("x value is:")
    disp(out.xList(:, out.hStop))
    disp("Relative errors are:")
    disp(out.relErr(:, out.hStop-1))
end

% convergence ratio: converge -> 16
convErr = zeros(length(x), size(out.approxErr, 2)-1);

for i = 1:size(convErr, 2)
    convErr(:, i) = out.approxErr(:, i)./out.approxErr(:, i+1);
end

% order estimate: p -> 4
p = log2(convErr);


% state vector / step size plots
figure("Name", "State Vector")
tiledlayout(2, 1)

nexttile
figA = loglog(1:nTests, out.xList, '--');
grid on
ylim padded
ax = gca;
ax.XTick = 0:1:nTests;
ylabel('state values')
hold on
lineFig = xline(out.hStop, 'r--', 'DisplayName', 'hStop', 'LineWidth', 1.2);
hold off

nexttile
plot(1:nTests, out.hList, 'x--', 'LineWidth', 1.2);
grid on
ylim padded
ax = gca;
ax.XTick = 0:1:nTests;
ylabel('step size')
xlabel('reduction step')
hold on
xline(out.hStop, 'r--', 'LineWidth', 1.2);
hold off

% generate legend
for i = 1:size(out.xList, 1)
    figA(i).DisplayName = "x" + i;
end
legend([figA; lineFig])

markers = ['o', '+', 'x', 's', 'd', '^', 'v', '*', '.'];
% generate larger markers array if required
markers = repmat(markers, [1, 1+round(length(x0)/length(markers))]);
for i = 1:length(figA)
    figA(i).Marker = markers(i);
    figA(i).LineWidth = 1.2;
end


% error / convergence plots
figure("Name", "Error and Convergence")
tiledlayout(4, 1)

nexttile(2)
title('rel error')
figures(2, :) = loglog(2:nTests, out.relErr);
grid on;
ylabel('rel error')
xlim([1 nTests+1])
ylim([-relTol/1000 relTol/100])
xl = xlim;
yl = ylim;
hold on
xline(out.hStop, 'r--', 'LineWidth', 1.2)
hold off

nexttile(1)
title('approx error')
figures(1, :) = loglog(2:nTests, out.approxErr);
grid on;
ylabel('approx error')
xlim(xl)
ylim(yl)
hold on
hStopFig = xline(out.hStop, 'r--', 'LineWidth', 1.2 , 'DisplayName', 'hStop');
hold off

nexttile
title ('convergence ratio')
figures(3, :) = plot(3:nTests, convErr);
ylabel('convergence ratio')
xlim(xl)
ylim([12 20])
hold on
xline(out.hStop, 'r--', 'LineWidth', 1.2)
yline(16, 'r')
hold off

nexttile
title('observed order')
figures(4, :) = plot(3:nTests, p);
xlabel('reduction step')
ylabel('observed order')
xlim(xl)
ylim([2 6])
hold on
xline(out.hStop, 'r--', 'LineWidth', 1.2)
yline(4, 'r')
hold off

% generate legend
for i = 1:size(out.approxErr, 1)
    figures(1, i).DisplayName = "x" + i;
end
nexttile(1)
legend([figures(1, :), hStopFig])

% generate marker and line styles
for i = 1:4
    fig = figures(i, :);
    
    for j = 1:length(fig)
        fig(j).Marker = markers(j);
    end

    set(fig, 'LineStyle', ':');
    set(fig, 'LineWidth', 1.5);
end

end