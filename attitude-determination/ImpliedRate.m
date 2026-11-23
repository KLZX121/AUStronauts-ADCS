function wImp = ImpliedRate(q, dt)
%ImpliedRate Body rate that exactly carries q(:, i-1) to q(:, i) over dt
%   wImp = ImpliedRate(q, dt)
%
%   Uses the same convention as MEKF.PropQDisc, i.e.
%       q(i) = q(i-1) (x) dq,   dq = [cos(|w|dt/2); sin(|w|dt/2) w/|w|]
%   with (x) the Hamilton product, so that propagating wImp(:, i) with
%   MEKF.PropQDisc reproduces the truth quaternion to machine precision.
%
%   Inputs
%   q       (4, n)  True attitude quaternion history (scalar first)
%   dt      (1, 1)  Timestep (s)
%
%   Output
%   wImp    (3, n)  Rate for the step ARRIVING at sample i (rad/s, body).
%                   Column 1 has no preceding step and is a copy of column 2.

n = size(q, 2);
wImp = zeros(3, n);

for i = 2:n
    qPrevConj = [q(1, i-1); -q(2:4, i-1)];
    dq = RIEKF.HamProd(qPrevConj, q(:, i));
    dq = dq/norm(dq);
    wImp(:, i) = RIEKF.QLog(dq)/dt;
end

wImp(:, 1) = wImp(:, 2);
end