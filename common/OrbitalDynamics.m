function xDot = OrbitalDynamics(x, d)
%OrbitalDynamics Returns the state ODEs for orbital r and v
%   xDot = OrbitalDynamics(x, d)
%
%   Inputs
%   x       (:, 1)  State vector
%   d       struct  Data struct
%   .iR     (1, 3)  Indices of orbital position in state vector (km)
%   .iV     (1, 3)  Indices of orbital velocity in state vector (km/s)
%   .mu     (1, 1)  Standard gravitational parameter (m^3/s^2)
%
%   Outputs
%   xDot    (6, 1)  [dr/dt; dv/dt] ODEs for r and v

r = x(d.iR);
v = x(d.iV);

mu = d.mu*1e-9;

rDot = v;
vDot = -mu/norm(r)^3 * r;

xDot = [rDot; vDot];

end