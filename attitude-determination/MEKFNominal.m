classdef MEKFNominal < MEKF
%MEKFNominal A 6-state MEKF for nominal mission
%   mekf = MEKFNominal
%
%   State Vector
%   xEst (6) = [
%           attitude error (3); 
%           gyro biases (3)
%   ]
%
%   See MEKF superclass for method documentation

properties(Constant)
    nStates = 6;
end
properties 
    C (:, :) double
end

methods
    % calculate the matrix C at initialisation
    function o = Initialise(o, x0, q0, P0, Q, gyro, mag, css)
        Initialise@MEKF(o, x0, q0, P0, Q, gyro, mag, css);

        o.C = eye(3) / (eye(3) + gyro.S);
    end
end
methods(Access=protected)
    function wEst = wEstFn(o, wGyro)
        bEst = o.xEst(4:6);

        wEst = o.C*(wGyro - bEst);
    end
    function [F, G] = StateMatrices(o)
        F = [
            -Skew(o.wEst), -o.C;
            zeros(3, 6);
            ];
        G = blkdiag(-o.C, eye(3));
    end
end
methods(Static)
    function Plot(t, xErrList, sigma3List)
        %Plot Plots the state errors over time, as well as 3-sigma bounds
        %   Plot(t, xErrList, sigma3List)
        %   
        %   Inputs
        %   t           (1, :)      List of timesteps (s)
        %   xErrList    (n, :)      List of error states
        %   sigma3List  (n, :)      List of 3-sigma bounds
        
        figData.figName = "NominalMEKF";
        figData.tlDim = [2, 1];
        figData.plotSizes = repmat([1 1], 2, 1);
        figData.yLabels = [
            "\delta\vartheta (deg)"; 
            "\Delta\beta (deg/h)"; 
        ];
        figData.y = {
            rad2deg(xErrList(1:3,:));
            rad2deg(xErrList(4:6, :)).*3600;
        };
        figData.sigma3s = {
            rad2deg(sigma3List(1:3, :));
            rad2deg(sigma3List(4:6, :)).*3600;
        };
        figData.xLabel = "t (min)";
        figData.x = t./60;

        MEKF.Plot(figData)
    end
end
end