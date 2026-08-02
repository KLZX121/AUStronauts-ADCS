function [bMag, hEst, HEst] = MagModel(bRef, q, mag, qEst, nStates)
%MagModel Magnetometer Model
%   bMag = MagModel(bRef, q, mag, qEst, nStates, M)
%
%   Simulates magnetometer measurements, and optionally returns measurement
%   matrices h and H for use in an MEKF when an estimated quaternion is
%   input. This is only applicable for MEKF's that do not estimate
%   magnetometer calibration parameters.
%   
%   Inputs
%   bRef            (3, 1)      true magnetic field - use IGRFECI.m (ECI)
%   q               (4, 1)      attitude quaternion
%   mag             (struct)    magnetometer calibration data
%   .sigma          (3, 1)      1-sigma noise for each axis
%   .bias           (3, 1)      biases
%   .D              (3, 3)      scale factor and non-orthogonality matrix
%   .O              (3, 3)      DCM for sensor to body frame rotation
%   .M              (3, 3)      matrix M = inv(eye(3) + D)
%   qEst            (4, 1)      estimated quaternion from MEKF (optional)
%   nStates         (1, 1)      number of states in MEKF (optional)
%   
%   Outputs
%   bMag            (3, 1)      vector magnetic field measurement (T) (body frame)
%   hEst            (3, 1)      estimated measurement function (optional)
%   HEst            (3, n)      estimated measurement matrix

A = QToDCM(q);

noise = mag.sigma.*randn(3, 1);

bMag = mag.M*(mag.O'*A*bRef + mag.bias) + noise;

% MEKF measurement functions and matrices
if (nargin > 3)
    mag.sigma = zeros(3, 1);
    bEst = MagModel(bRef, qEst, mag);

    hEst = bEst;

    HEst = [
        mag.M*mag.O'*Skew(QToDCM(qEst)*bRef), ...
        zeros(3, nStates-3)
    ];
end

end
