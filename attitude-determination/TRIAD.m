function [A, P] = TRIAD(b, r, sigmaMag)
%TRIAD Implement the TRIAD attitude determination algorithm
%   [A, P] = TRIAD(b, r)
%
%   Transforms reference -> body
%   (src - Wertz 1978)
%
%   Inputs
%   b           (3, 2)  [b1, b2] Body frame unit vectors (sensor measurements)
%   r           (3, 2)  [r1, r2] Reference frame unit vectors
%   sigmaMag    (1, 1)  Magnetometer 1-sigma noise (T)
%
%   Outputs
%   A           (3, 3)  Attitude DCM
%   P           (3, 3)  Covariance matrix

b1 = b(:, 1);
b2 = b(:, 2);
r1 = r(:, 1);
r2 = r(:, 2);

% body triad
w1 = b1;

w2 = cross(b1, b2);
w2 = w2/norm(w2);

w3 = cross(w1, w2);

M_b = [w1 w2 w3];

% reference triad
v1 = r1;

v2 = cross(r1, r2);
v2 = v2/norm(v2);

v3 = cross(v1, v2);

M_r = [v1 v2 v3];

% DCM ref -> body
A = M_b*M_r';

% covariance
sigma1 = sigmaMag/norm(b1);
% sigma2 determined empirically through simulation
sigma2 = 0.07;
P = sigma1^2.*eye(3) ...
        + (1/norm(cross(b1, b2))^2).*((sigma2^2 - sigma1^2).*b1.*b1' ...
        + sigma1^2*dot(b1, b2).*(b1*b2' + b2*b1'));

end