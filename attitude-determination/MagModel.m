classdef MagModel < handle
properties
    sigma (3, 1)
    bias (3, 1)
    D (3, 3)
    O (3, 3)

    M (3, 3)
end

methods
    function o = MagModel(bias, D, O, sigma)
        %MagModel Models a 3-axis magnetometer
        %   o = MagModel(bias, D, O, sigma) 
        %   
        %   Inputs
        %   bias    (3, 1)      biases (T)
        %   D       (3, 3)      scale factor and non-orthogonality matrix
        %   O       (3, 3)      DCM for sensor to body frame rotation
        %   sigma   (3, 1)      1-sigma noise for each axis (T)
        %   
        %   Outputs
        %   o       (object)    MagModel object       

        o.bias = bias;
        o.D = D;
        o.O = O;
        o.sigma = sigma;
        
        o.M = eye(3) / (eye(3) + D);
    end

    function bMag = Measurement(o, bRef, q, genNoise)
        %Measurement Simulates magnetometer measurements
        %   bMag = Measurement(o, bRef, q, genNoise)
        %
        %   Inputs
        %   bRef        (3, 1)  True magnetic field in reference frame (T)
        %   q           (4, 1)  True attitude quaternion
        %   genNoise    (bool)  Whether to generate noise
        %
        %   Outputs
        %   bMag        (3, 1)  Measured magnetic field in body frame (T)

        A = QToDCM(q);

        if (genNoise)
            noise = o.sigma.*randn(3, 1);
        else
            noise = zeros(3, 1);
        end

        bMag = o.M*(o.O'*A*bRef + o.bias) + noise;
    end

    function [h, H] = MEKFMatrices(o, bRef, qEst, nStates)
        %MEKFMatrices Returns functions and matrices for a MEKF
        %   [h, H] = MEKFMatrices(o, bRef, qEst, nStates)
        %   
        %   Returns the estimated measurement function h and measurement
        %   sensitivity matrix H, assuming the MEKF does not estimate
        %   magnetometer calibration parameters
        %
        %   Inputs
        %   bRef        (3, 1)      Reference magnetic field (ECI)
        %   qEst        (4, 1)      Estimated attitude quaternion
        %   nStates     (double)    Number of states in filter
        %
        %   Outputs
        %   h           (3, 1)      Estimated measurement (T)
        %   H           (3, n)      Measurement sensitivity matrix

        h = o.Measurement(bRef, qEst, false);

        H = [
            o.M*o.O'*Skew(QToDCM(qEst)*bRef), ...
            zeros(3, nStates-3)
            ];
    end
end
end