function [x, xList] = PropState(xDot, x0, d, h, t)
%PropState Propagates a state vector with RK4
%   [x, xList] = PropState(xDot, x0, d, h, t0, tf)
%   
%   Inputs
%   xDot    (n, 1)      State variable ODE function handle
%   x0      (n, 1)      Initial state vector
%   d       (struct)    Data structure of parameters to pass to ODE
%   h                   Step size
%   t       (:, 1)      Time steps for propagation
%   
%   Outputs
%   x       (n, 1)      Final state vector
%   xList   (n, :)      List of state vectors at each timestep

N = (t(end)-t(1))/h;

xList = zeros(size(x0, 1), N+1);
xList(:, 1) = x0;

x = x0;
for i = 1:N
    x = RK4(xDot, x, h, t(i), d);
    xList(:, i+1) = x;
end

end