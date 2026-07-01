clear;
close all;

% specific fns, values
x0 = [0; 0; 0;];
xDotFn = @(x, t, T, I) [EulerDynamics(x(1:3), T, I)];
I = InertiaCubeSat('3U', 6);
T = [1; 1; 1;];




x = x0;

% test time interval
t0 = 0;
tf = 1;

% number of tests
n_tests = 17;

x_list = zeros(n_tests, length(x));
h_list = zeros(n_tests, 1);

% initial step
h_list(1) = 1;

% reduction formula
h_new = @(h) 0.5*h;

for test = 1:n_tests
    h = h_list(test);

    t = t0:h:tf;
    N = (tf-t0)/h;
    
    x = x0;
    
    L = zeros(N, 3);
    for i = 1:N
        x = RK4(xDotFn, x, h, t(i), T, I);

        L(i, :) = (I*x)';
    end
    x_list(test, :) = x;

    h_list(test+1) = h_new(h);
end

x_list
h_list




% absolute error: converge -> 0
err = zeros(n_tests-1, length(x));

for i = 1:n_tests-1
    err(i, :) = abs(x_list(i+1, :) - x_list(i, :));
end

err

title('abs error')
plot(2:n_tests, err, '*--');
grid on;
xlabel('reduction step')
ylabel('abs error')
ylim([0 1])

hold on
yline(0, 'r')
hold off

% convergence ratio: converge -> 16
conv_err = zeros(length(err)-1, length(x));

for i = 1:length(conv_err)
    conv_err(i, :) = err(i, :)./err(i+1, :);
end

conv_err

figure;
title ('convergence ratio')
plot(3:n_tests, conv_err, '*--');
grid on;
xlabel('reduction step')
ylabel('convergence ratio')
ylim([12 20])

hold on
yline(16, 'r')
hold off

% order estimate: p -> 4
p = log(conv_err)/log(2)

figure;
title('observed order')
plot(3:n_tests, p, '*--')
grid on;
xlabel('reduction step')
ylabel('observed order')
ylim([2 6])

hold on
yline(4, 'r')
hold off