classdef MEKFGyroCal < MEKF
%MEKFGyroCal A 15-state gyro calibration MEKF
%   mekf = MEKFGyroCal
%
%   State Vector
%   xEst (15) = [
%           attitude error (3); 
%           bias (3); 
%           scale factors (3); 
%           misalignment upper (3); 
%           misalignment lower (3)
%   ]
%
%   See MEKF superclass for method documentation

properties(Constant)
    nStates = 15;
end
properties
    CEst (:, :) double
end

methods(Access=protected)
    function wEst = wEstFn(o, wGyro)
        bEst = o.xEst(4:6);

        sEst = o.xEst(7:9);
        kUEst = o.xEst(10:12);
        kLEst = o.xEst(13:15);
        SEst = GyroModel.SMatrix(sEst, kUEst, kLEst);
        o.CEst = eye(3) / (eye(3) + SEst);

        wEst = o.CEst*(wGyro - bEst);
    end
    function [F, G] = StateMatrices(o)
        U = [
            o.wEst(2), o.wEst(3), 0;
            0, 0, o.wEst(3);
            0, 0, 0
            ];
        L = [
            0, 0, 0;
            o.wEst(1), 0, 0;
            0, o.wEst(1), o.wEst(2);
            ];

        F = [
            -Skew(o.wEst), -o.CEst, -o.CEst*diag(o.wEst), -o.CEst*U, -o.CEst*L;
            zeros(12, 15)
            ];
        G = blkdiag(-o.CEst,eye(12));
    end
end
methods(Static)
    function Plot(t, xErrList, sigma3List, plotType)
        %Plot Plots the state errors over time, as well as 3-sigma bounds
        %   Plot(t, xErrList, sigma3List, plotType)
        %   
        %   Inputs
        %   t           (1, :)      List of timesteps (min)
        %   xErrList    (n, :)      List of error states
        %   sigma3List  (n, :)      List of 3-sigma bounds
        %   plotType    (str)       'overview' or 'detailed'

        yLabels = [
            "\delta\vartheta (deg)";
            "\Delta\beta (deg/h)";
            "\Deltas"; 
            "\Deltak_U"; 
            "\Deltak_L"
        ];
        y = {
            rad2deg(xErrList(1:3,:));
            rad2deg(xErrList(4:6, :)).*3600;
            xErrList(7:9, :);
            xErrList(10:12, :);
            xErrList(13:15, :);
        };
        sig3 = {
            rad2deg(sigma3List(1:3, :));
            rad2deg(sigma3List(4:6, :)).*3600;
            sigma3List(7:9, :);
            sigma3List(10:12, :);
            sigma3List(13:15, :);
        };

        if plotType == "overview"
            figData.figName = "MEKF (15-state)";
            figData.plotType = plotType;

            figData.tlDim = [3, 2];
            figData.plotSizes = [[1 2]; repmat([1 1], 4, 1)];
            figData.yLabels = yLabels;
            figData.y = y;
            figData.sigma3s = sig3;

            figData.xLabel = "t (min)";
            figData.x = t;

            MEKF.Plot(figData)
        elseif plotType == "detailed"
            figNames = [
                "MEKF Attitude Error";
                "MEKF Gyro Bias Error";
                "MEKF Gyro Scale Factor Error";
                "MEKF Gyro Upper Misalignment Error";
                "MEKF Gyro Lower Misalignment Error";
            ];

            for i = 1:length(figNames)
                figData.figName = figNames(i);
                figData.plotType = plotType;
    
                figData.tlDim = [3, 1];
                figData.plotSizes = repmat([1 1], 3, 1);
                figData.yLabels = [
                    sprintf("X %s", yLabels(i));
                    sprintf("Y %s", yLabels(i))
                    sprintf("Z %s", yLabels(i))
                ];
                figData.y = {
                    y{i}(1, :);
                    y{i}(2, :);
                    y{i}(3, :);
                };
                figData.sigma3s = {
                    sig3{i}(1, :);
                    sig3{i}(2, :);
                    sig3{i}(3, :);
                };
    
                figData.xLabel = "t (min)";
                figData.x = t;
    
                MEKF.Plot(figData)
            end
        end
    end
end
end