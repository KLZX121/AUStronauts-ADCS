classdef CSSModel < handle
properties
    sigmaTheta (1, 1) double
    sigmaDark (1, 1) double
    sigmaEdge (1, 1) double
    sigmaYFn

    nSensors (1, 1) double
    uSensors (3, :) double
    fov (1, 1) double
    thetaFov (1, 1) double
    cosFov (1, 1) double
    yFov (1, 1) double

    LUT (:, 2) double
    yFn (1, 1) struct
    thetaFn (1, 1) struct
    dydthetaFn (1, 1) struct

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
        %   .sigmaTheta     (1, 1)      1-sigma angular measurement noise inside 
        %                               sensor FOV (rad)
        %   .sigmaEdge      (1, 1)      1-sigma angular noise outside sensor FOV 
        %                               (rad)
        %   .sigmaDark      (1, 1)      1-sigma angular noise for unilluminated
        %                               sensor (rad)
        %   .yLims          (1, 2)      [yMin, yMax] which is the minimum and
        %                               maximum possible measurement value (V or I)
        %   .T0             (1, 1)      Designated temperature of LUT (degC)
        %   .alpha          (1, 1)      Temperature coefficient (V or I) / degC
        %
        %   Outputs
        %   o               (obj)       CSSModel object

        o.nSensors = calCSS.nSensors;
        o.uSensors = calCSS.uSensors;

        o.yMin = calCSS.yLims(1);
        o.yMax = calCSS.yLims(2);
        o.T0 = calCSS.T0;
        o.alpha = calCSS.alpha;

        o.fov = calCSS.fov;
        o.thetaFov = o.fov/2;
        o.cosFov = cos(o.thetaFov);

        o.sigmaTheta = calCSS.sigmaTheta;
        o.sigmaEdge = calCSS.sigmaEdge;
        o.sigmaDark = calCSS.sigmaDark;

        % create LUT interpolation functions
        o.LUT = calCSS.LUT;
        % y = f(theta) (piecewise polynomial struct)
        o.yFn = pchip(o.LUT(:, 1), o.LUT(:, 2));
        % theta = f(y)
        o.thetaFn = pchip(o.LUT(:, 2), o.LUT(:, 1));

        o.yFov = ppval(o.yFn, o.thetaFov);

        % measurement noise calculation
        grad = gradient(o.LUT(:, 2), o.LUT(:, 1));
        % y' = f'(theta)
        o.dydthetaFn = pchip(o.LUT(:, 1), grad);
        % sigmaY function for fov measurements
        o.sigmaYFn = @(theta) abs(ppval(o.dydthetaFn, theta)).*o.sigmaTheta;
    end

    function [y, yLit] = Measurement(o, q, uSunRef, T, genNoise)
        %Measurement Simulates sun sensor measurements
        %   [y, yLit] = Measurement(o, q, uSunRef, T, genNoise)
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
        %   yLit        (n, 1)  Indicates whether sensor is:
        %                       lit in fov = 1
        %                       lit outside fov = 2
        %                       unlit = 0

        % rotate sun reference vector to body frame
        uSunBody = QToDCM(q)*uSunRef;
        
        % get true incidence angle      
        cosThetas = clip(o.uSensors' * uSunBody, -1, 1);
        trueThetas = acos(cosThetas);

        % temp offset
        yTOffset = o.alpha*(T - o.T0);

        % sensor illumination status
        yLit = zeros(6, 1);

        % get measurement including temp offset
        % TODO: add solar panel occlusion (and fgm??)
        y = zeros(o.nSensors, 1);
        if (genNoise) 
            sigmaY = zeros(o.nSensors, 1);
        end
        for i = 1:o.nSensors
            if (cosThetas(i) >= o.cosFov)
                % if sun is within fov

                % use LUT and pchip interpolation
                y(i) = ppval(o.yFn, trueThetas(i));

                % add temperature offset
                y(i) = y(i) + yTOffset;

                % measurement noise (using datasheet angular error)
                if (genNoise)
                    sigmaY(i) = o.sigmaYFn(trueThetas(i));
                end

                yLit(i) = 1;
            elseif (cosThetas(i) > 0)
                % if sun is outside fov but within 90 deg

                % find amplitude of cosine by generating continuous
                % curve from the fov limit of the LUT
                amp = (o.yFov + yTOffset) / o.cosFov;

                y(i) = amp*cosThetas(i);

                % measurement noise (higher than LUT)
                if (genNoise)
                    sigmaY(i) = amp*sin(trueThetas(i))*o.sigmaEdge;
                end

                yLit(i) = 2;
            else
                % sun is towards back-side of sensor
                y(i) = 0;

                if (genNoise)
                    sigmaY(i) = o.sigmaDark;
                end

                yLit(i) = 0;
            end
        end

        % add noise
        if (genNoise)
            y = y + sigmaY.*randn(o.nSensors, 1);
        end

        % clip measurement to physical y limits
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

        yTOffset = o.alpha*(T - o.T0);

        yFovMeas = o.yFov + yTOffset;

        yThreshold = 3*o.sigmaDark;

        c = zeros(o.nSensors, 1);
        for i = 1:o.nSensors
            if (y(i) <= yThreshold)
                c(i) = 0;

            elseif (y(i) > yFovMeas)
                % measurement in fov
                % use lut to get angle

                yLUT = y(i) - yTOffset;
                yLUT = clip(yLUT, min(o.LUT(:, 2)), max(o.LUT(:, 2)));

                t = ppval(o.thetaFn, yLUT);

                c(i) = cos(t);
            elseif (y(i) > 0)
                % measurement outside fov
                % use reverse cos
                
                c(i) = y(i)/(yFovMeas/o.cosFov);
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