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
%           .hI         Initial step (default 1)
%
%   Outputs:
%   out     (struct)    Ouput struct
%           .hStop                      Step index at which relative error 
%                                       is below relative tolerance
%                                       (nSigFig)
%           .xList      (nTests, nx)    List of states for each test
%           .hList      (nTests, 1)     List of steps
%           .approxErr  (nTests-1, nx)  Approximate error (x_i - x_(i-1))
%           .relErr     (nTests-1, nx)  Relative error (approxErr/x_i)

if (nargin < 5)
    nSigFig = 3;
end
if (nargin < 6)
    dTest.t0 = 0;
    dTest.tf = 1;
    dTest.hI = 1;
end

% good implementation of rk4 in chapra pg. 638
% verify with ode45

out.xList = zeros(nTests, length(x0));
out.hList = zeros(nTests, 1);
out.hList(1) = dTest.hI;

% reduction formula
hNew = @(h) 0.5*h;

for testI = 1:nTests
    % get step
    h = out.hList(testI);

    % run RK4 with step
    t = dTest.t0:h:dTest.tf;
    N = (dTest.tf-dTest.t0)/h;
    
    x = x0;
    
    for i = 1:N
        x = RK4(xDotFn, x, h, t(i), d);
    end
    out.xList(testI, :) = x;

    % generate next step if not last step already
    if (testI < nTests)
        out.hList(testI+1) = hNew(h);
    end
end

% approx error, relative error: converge -> 0
out.approxErr = zeros(nTests-1, length(x));
out.relErr = zeros(size(out.approxErr));

% stopping criterion
% note relTol is in % while relErr is not.
% src - chapra pg. 104
relTolFn = @(n) 0.5*10^(2-n);
relTol = relTolFn(nSigFig);
disp('Least amount of accurate sig figs is ' + string(nSigFig))
disp('Relative tolerance is ' + string(relTol/100));

% step at which relTol is met
out.hStop = 0;

for i = 1:nTests-1
    out.approxErr(i, :) = abs(out.xList(i+1, :) - out.xList(i, :));
    out.relErr(i, :) = abs(out.approxErr(i, :) ./ out.xList(i+1, :));

    if ((out.relErr(i, :)*100 < relTol) & (out.hStop == 0))  
        out.hStop = i+1;
        disp('Tolerance met at step ' + string(out.hStop));
        disp('Step size is ' + string(out.hList(out.hStop)))
        disp("x value is " + mat2str(out.xList(out.hStop, :)))
    end
end
if (~out.hStop)
    disp("No appropriate step was found - manually check values and plots")
end

% convergence ratio: converge -> 16
convErr = zeros(length(out.approxErr)-1, length(x));

for i = 1:length(convErr)
    convErr(i, :) = out.approxErr(i, :)./out.approxErr(i+1, :);
end

% order estimate: p -> 4
p = log2(convErr);


% state vector / step size plots
figure("Name", "State Vector")
tiledlayout(2, 1)

nexttile
figA = plot(1:nTests, out.xList, '--');
grid on
ylim padded
ax = gca;
ax.XTick = 0:1:nTests;
ylabel('state values')
hold on
lineFig = xline(out.hStop, 'k', 'DisplayName', 'hStop');
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
xline(out.hStop, 'k');
hold off

% generate legend
for i = 1:size(out.xList, 2)
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
figures(2, :) = plot(2:nTests, out.relErr);
grid on;
ylabel('rel error')
ylim([-relTol/1000 relTol/100])
xl = xlim;
yl = ylim;
hold on
xline(out.hStop, 'k')
hold off

nexttile(1)
title('approx error')
figures(1, :) = plot(2:nTests, out.approxErr);
grid on;
ylabel('approx error')
ylim(yl)
hold on
hStopFig = xline(out.hStop, 'k', 'DisplayName', 'hStop');
hold off

nexttile
title ('convergence ratio')
figures(3, :) = plot(3:nTests, convErr);
grid on;
ylabel('convergence ratio')
xlim(xl)
ylim([12 20])
hold on
xline(out.hStop, 'k')
yline(16, 'r')
hold off

nexttile
title('observed order')
figures(4, :) = plot(3:nTests, p);
grid on;
xlabel('reduction step')
ylabel('observed order')
xlim(xl)
ylim([2 6])
hold on
xline(out.hStop, 'k')
yline(4, 'r')
hold off

% generate legend
for i = 1:size(out.approxErr, 2)
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