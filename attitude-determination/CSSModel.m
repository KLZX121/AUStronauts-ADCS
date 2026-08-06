classdef CSSModel < handle
properties
    sigma (1, 1) double
    sigmaTheta (1, 1) double

    nSensors (1, 1) double
    uSensors (3, :) double
    fov (1, 1) double

    LUT (:, 2) double
    yFn (1, 1) struct
    thetaFn (1, 1) struct

    yMin (1, 1) double
    yMax (1, 1) double
    T0 (1, 1) double
    alpha (1, 1) double
end

methods
    function o = CSSModel(calCSS)
        %CSSModel Simulates a set of coarse sun sensors
        %   o = CSSModel(calCSS)
        %
        %   Inputs
        %   calCSS          (struct)    Calibration Parameters
        %   .nSensors       (1, 1)      Number of sun sensors
        %   .uSensors       (3, n)      Unit normal vectors of each sensor
        %   .fov            (1, 1)      Total field of view of sensors (rad)
        %   .LUT            (:, 2)      Look Up Table of angles (rad) in column 1
        %                               and sensor output (V or I) in column 2
        %   .sigma          (1, 1)      1-sigma measurement noise of sensors (V or I)
        %   .sigmaTheta     (1, 1)      1-sigma angular measurement noise (rad)
        %   .yLims          (1, 2)      [yMin, yMax] which is the minimum and
        %                               maximum possible measurement value (V or I)
        %   .T0             (1, 1)      Designated temperature of LUT (degC)
        %   .alpha          (1, 1)      Temperature coefficient (V or I) / degC
        %
        %   Outputs
        %   o               (obj)       CSSModel object

        o.nSensors = calCSS.nSensors;
        o.uSensors = calCSS.uSensors;
        o.fov = calCSS.fov;
        o.sigma = calCSS.sigma;
        o.sigmaTheta = calCSS.sigmaTheta;
        o.yMin = calCSS.yLims(1);
        o.yMax = calCSS.yLims(2);
        o.T0 = calCSS.T0;
        o.alpha = calCSS.alpha;

        % create interpolation functions
        o.LUT = calCSS.LUT;
        % y = f(theta) (piecewise polynomial struct)
        o.yFn = pchip(calCSS.LUT(:, 1), calCSS.LUT(:, 2));
        % theta = f(y)
        o.thetaFn = pchip(calCSS.LUT(:, 2), calCSS.LUT(:, 1));
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
                y(i) = ppval(o.yFn, trueThetas(i));

                % add temperature offset
                y(i) = y(i) + yTOffset;
            elseif (cosThetas(i) > 0)
                % if sun is outside fov but within 90 deg

                % find amplitude of cosine by generating continuous
                % curve from the fov limit of the LUT
                yFov = ppval(o.yFn, thetaFov);
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

        yFov = ppval(o.yFn, thetaFov);
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

                t = ppval(o.thetaFn, yLUT);

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

    function [h, H] = MEKFMatrices(o, qEst, uSunRef, T, nStates, sensorI)
        %MEKFMatrices Calculates estimated measurements for MEKF
        %   [h, H] = MEKFMatrices(o, qEst, uSunRef, T, nStates, sensorI)
        %
        %   Returns matrices for all sensors or one sensor if sensorI is
        %   specified.
        %
        %   Inputs
        %   qEst        (4, 1)  Estimated quaternion
        %   uSunRef     (3, 1)  Unit true reference Sun vector (ECI)
        %   T           (1, 1)  Current temperature (degC)
        %   nStates     (1, 1)  Number of filter states
        %   sensorI     (1, 1)  Index of sensor to get (optional)
        %
        %   Outputs
        %   h           (:, 1)  Estimated measurement
        %   H           (:, :)  Measurement sensitivity matrix

        h = o.Measurement(qEst, uSunRef, T, false);

        uCSS = o.uSensors;
        
        if (nargin > 5)
            uCSS = uCSS(:, sensorI);
            h = h(sensorI);
        end
        
        uSunEst = QToDCM(qEst)*uSunRef;
        
        thetaEst = acos(uCSS' * uSunEst);

        dTheta = deg2rad(0.01);
        yPlus = ppval(o.yFn, thetaEst + dTheta);
        yMinus = ppval(o.yFn, thetaEst - dTheta);

        dy = (yPlus - yMinus)/(2*dTheta);

        H = -dy.*(uCSS'*Skew(uSunEst))./sin(thetaEst);

        if (nargin > 5)
            H = [H zeros(1, nStates-3)];
        else
            H = [H zeros(o.nSensors, nStates-3)];
        end
    end
end
end