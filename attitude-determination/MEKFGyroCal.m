classdef MEKFGyroCal < MEKF
%MEKFGyroCal A 15-state gyro calibration MEKF
%   mekf = MEKFGyroCal
%
%   State Vector
%   x(15) = [
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

methods(Static)
    function Plot(t, xErrList, sigma3List)
        %Plot Plots the state errors over time, as well as 3-sigma bounds
        %   Plot(t, xErrList, sigma3List)
        %   
        %   Inputs
        %   t           (1, :)      List of timesteps (s)
        %   xErrList    (n, :)      List of error states
        %   sigma3List  (n, :)      List of 3-sigma bounds

        figData.figName = "Gyro Calibration MEKF";
        figData.tlDim = [3, 2];
        figData.plotSizes = [[1 2]; repmat([1 1], 4, 1)];
        figData.yLabels = [
            "\delta\vartheta (deg)"; 
            "\Delta\beta (deg/h)"; 
            "\Deltas"; 
            "\Deltak_U"; 
            "\Deltak_L"
        ];
        figData.y = {
            rad2deg(xErrList(1:3,:));
            rad2deg(xErrList(4:6, :)).*3600;
            xErrList(7:9, :);
            xErrList(10:12, :);
            xErrList(13:15, :);
        };
        figData.sigma3s = {
            rad2deg(sigma3List(1:3, :));
            rad2deg(sigma3List(4:6, :)).*3600;
            sigma3List(7:9, :);
            sigma3List(10:12, :);
            sigma3List(13:15, :);
        };
        figData.xLabel = "t (min)";
        figData.x = t./60;

        MEKF.Plot(figData)
    end
end
methods(Access=protected)
    function wEst = wEstFn(o, wGyro)
        bEst = o.xEst(4:6);

        sEst = o.xEst(7:9);
        kUEst = o.xEst(10:12);
        kLEst = o.xEst(13:15);
        SEst = MEKF.SMatrix(sEst, kUEst, kLEst);
        o.CEst = eye(3) / (eye(3) + SEst);

        wEst = o.CEst*(wGyro - bEst);
    end
    function [y, h, H] = MeasurementMatrices(o, extData)
        noise = extData.sigmaST*randn(3, 1);
        dq = [1; 0.5*noise];
        dq = dq./norm(dq);
        qStar = QProd(dq, extData.qTrue);
    
        qMeas = QProd(qStar, QConj(o.qEst));
        y = 2*qMeas(2:4)./qMeas(1);
    
        % measurement sensitivity matrix
        h = zeros(3, 1);
        H = [eye(3) zeros(3, o.nStates-3)];
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
end