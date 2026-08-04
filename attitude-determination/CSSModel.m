classdef CSSModel < handle
properties
    sigma (1, 1) double
    sigmaTheta (1, 1) double

    nSensors (1, 1) double
    uSensors (3, :) double
    fov (1, 1) double
    LUT (:, 2) double

    yMin (1, 1) double
    yMax (1, 1) double
    T0 (1, 1) double
    alpha (1, 1) double
end

methods
    function o = CSSModel(nSensors, uSensors, fov, LUT, sigma, sigmaTheta, yLims, T0, alpha)
        %CSSModel Simulates a set of coarse sun sensors
        %   o = CSSModel(nSensors, uSensors, fov, LUT, sigma, yLims, T0, alpha)
        %
        %   Inputs
        %   nSensors    (1, 1)  Number of sun sensors
        %   uSensors    (3, n)  Unit normal vectors of each sensor
        %   fov         (1, 1)  Total field of view of sensors (rad)
        %   LUT         (:, 2)  Look Up Table of angles (rad) in column 1
        %                       and sensor output (V or I) in column 2
        %   sigma       (1, 1)  1-sigma measurement noise of sensors (V or I)
        %   sigmaTheta  (1, 1)  1-sigma angular measurement noise (rad)
        %   yLims       (1, 2)  [yMin, yMax] which is the minimum and
        %                       maximum possible measurement value (V or I)
        %   T0          (1, 1)  Designated temperature of LUT (degC)
        %   alpha       (1, 1)  Temperature coefficient (V or I) / degC
        %
        %   Outputs
        %   o           (obj)   CSSModel object

        o.nSensors = nSensors;
        o.uSensors = uSensors;
        o.fov = fov;
        o.LUT = LUT;
        o.sigma = sigma;
        o.sigmaTheta = sigmaTheta;
        o.yMin = yLims(1);
        o.yMax = yLims(2);
        o.T0 = T0;
        o.alpha = alpha;
    end

    function y = Measurement(o, q, uSunRef, T, genNoise)
        %Measurement Simulates sun sensor measurements
        %   y = Measurement(o, q, uSunRef, T, genNoise)
        %
        %   Simulates measurements considering temperature and noise using 
        %   a LUT. Does not consider interference sources like albedo
        %   
        %   Inputs
        %   q           (4, 1)  True attitude quaternion
        %   uSunRef     (3, 1)  True unit sun vector in reference frame (ECI)
        %   T           (1, 1)  Current temperature (degC)
        %   genNoise    (bool)  Whether to generate measurement noise
        %
        %   Outputs
        %   y           (n, 1)  Measurements of each sun sensor (V or I)

        % rotate sun reference vector to body frame
        uSunBody = QToDCM(q)*uSunRef;
        
        cosThetas = clip(o.uSensors' * uSunBody, -1, 1);
        trueThetas = acos(cosThetas);

        % add angular noise
        if (genNoise)
            thetaNoise = o.sigmaTheta.*randn(o.nSensors, 1);
            trueThetas = trueThetas + thetaNoise;

            cosThetas = clip(cos(trueThetas), -1, 1);
        end

        thetaFov = o.fov/2;
        cosFov = cos(thetaFov);

        % temp offset
        yTOffset = o.alpha*(T - o.T0);

        % ideal measurement
        % TODO: add solar panel occlusion (and fgm??)
        y = zeros(o.nSensors, 1);
        for i = 1:o.nSensors
            if (cosThetas(i) >= cosFov)
                % if sun is within fov

                % use LUT and pchip interpolation
                y(i) = interp1(o.LUT(:, 1), o.LUT(:, 2), trueThetas(i), 'pchip');

                % add temperature offset
                y(i) = y(i) + yTOffset;
            elseif (cosThetas(i) > 0)
                % if sun is outside fov but within 90 deg

                % find amplitude of cosine by generating continuous
                % curve from the fov limit of the LUT
                yFov = interp1(o.LUT(:, 1), o.LUT(:, 2), thetaFov, 'pchip');
                amp = (yFov + yTOffset) / cosFov;

                y(i) = amp*cosThetas(i);
            else
                % sun is towards back-side of sensor
                y(i) = 0;
            end
        end

        % add measurement noise
        if (genNoise)
            y = y + o.sigma.*randn(o.nSensors, 1);
        end

        % clip measurement
        y = clip(y, o.yMin, o.yMax);
    end

    function uS = CalcSunVec(o, y, T)
        %CalcSunVec Calculates the sun vector from 6 orthogonal CSS
        %   uS = CalcSunVec(o, y, T)
        %
        %   Assumes o.uSensors = [+x; +y; +z; -x; -y; -z];
        %
        %   Inputs
        %   y   (6, 1)  Sun sensor measurements (from o.Measurement)
        %   T   (1, 1)  Current temperature (degC)
        %
        %   Outputs
        %   uS  (3, 1)  Unit Sun vector in body frame

        thetaFov = o.fov/2;
        cosFov = cos(thetaFov);

        yTOffset = o.alpha*(T - o.T0);

        yFov = interp1(o.LUT(:, 1), o.LUT(:, 2), thetaFov, 'pchip');
        yFov = yFov + yTOffset;

        yThreshold = 3*o.sigma;

        c = zeros(o.nSensors, 1);
        for i = 1:o.nSensors
            if (y(i) <= yThreshold)
                c(i) = 0;

            elseif (y(i) > yFov)
                % measurement in fov
                % use lut to get angle

                yLUT = y(i) - yTOffset;
                yLUT = clip(yLUT, min(o.LUT(:, 2)), max(o.LUT(:, 2)));

                t = interp1(o.LUT(:, 2), o.LUT(:, 1), yLUT, 'pchip');

                c(i) = cos(t);
            elseif (y(i) > 0)
                % measurement outside fov
                % use reverse cos
                
                c(i) = y(i)/(yFov/cosFov);
            end
        end
        c = clip(c, 0, 1);

        uS = c(1:3) - c(4:6);
        uS = uS./norm(uS);
    end

    function [h, H] = MEKFMatrices(o, qEst, uSunRef, T, nStates)
        h = o.Measurement(qEst, uSunRef, T, false);

        H = 0;

    end
end
end