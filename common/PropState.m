function x = PropState(xDot, x0, d, h)
%PropState Propagates a state vector with RK4 by one timestep for an
%          autonomous system
%   x = PropState(xDot, x0, d, h)
%
%   Propagates a time-invariant (autonomous) system by one RK4 step
%   
%   Inputs
%   xDot    (n, 1)      State variable ODE function handle @(x, d)
%   x0      (n, 1)      Initial state vector
%   d       (struct)    Data structure of parameters to pass to ODE
%   h                   Step size for time
%   
%   Outputs
%   x       (n, 1)      Final state vector

h2 = h/2;

k1 = xDot(x0, d);
k2 = xDot(x0 + k1*h2, d);
k3 = xDot(x0 + k2*h2, d);
k4 = xDot(x0 + k3*h, d);

x = x0 + (h/6)*(k1 + 2*(k2 + k3) + k4);

end