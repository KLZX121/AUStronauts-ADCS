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
methods(Static)
    function Plot(t, xErrList, sigma3List)
        function bound = yBounds(sigVals, xVals)
            % sets smart y limits
            
             % get the timestep with smallest sigma-3 bound
            [~, minI] = min(sum(sigVals));

            % get max of [sigma-3; x] at this timestep
            % set y limits to 3 times this value
            y = 3*(max([sigVals(:, minI); xVals(:, minI)]));

            bound = [-y y];
        end

        figure;
        tl = tiledlayout(3, 2);
        tl.Title.String = "MEKF State Errors";
        colororder(lines(3))
        set(0, 'DefaultLineLineWidth', 1.4)
        
        nexttile([1 2])
        vals = rad2deg(xErrList(1:3,:));
        sigmaVals = rad2deg(sigma3List(1:3, :));
        plot(t./60, vals)
        grid on
        ylabel('\delta\vartheta (deg)')
        xlabel('t (min)')
        hold on
        plot(t./60, sigmaVals, ':')
        plot(t./60, -sigmaVals, ':')
        hold off
        ylim(yBounds(sigmaVals, vals));
        
        lg = legend('x', 'y', 'z', '3-\sigma_x', '3-\sigma_y', '3-\sigma_z');
        lg.Location = "eastoutside";
        
        nexttile
        vals = rad2deg(xErrList(4:6, :)).*3600;
        sigmaVals = rad2deg(sigma3List(4:6, :)).*3600;
        plot(t./60, vals)
        grid on
        ylabel('\Delta\beta (deg/h)')
        xlabel('t (min)')
        hold on
        plot(t./60, sigmaVals, ':')
        plot(t./60, -sigmaVals, ':')
        hold off
        ylim(yBounds(sigmaVals, vals));
        
        nexttile
        vals = xErrList(7:9, :);
        sigmaVals = sigma3List(7:9, :);
        plot(t./60, vals)
        grid on
        ylabel('\Deltas')
        xlabel('t (min)')
        hold on
        plot(t./60, sigmaVals, ':')
        plot(t./60, -sigmaVals, ':')
        hold off
        ylim(yBounds(sigmaVals, vals));
        
        nexttile
        vals = xErrList(10:12, :);
        sigmaVals = sigma3List(10:12, :);
        plot(t./60, vals)
        grid on
        ylabel('\Deltak_U')
        xlabel('t (min)')
        hold on
        plot(t./60, sigmaVals, ':')
        plot(t./60, -sigmaVals, ':')
        hold off
        ylim(yBounds(sigmaVals, vals));
        
        nexttile
        vals = xErrList(13:15, :);
        sigmaVals = sigma3List(13:15, :);
        plot(t./60, vals)
        grid on
        ylabel('\Deltak_L')
        xlabel('t (min)')
        hold on
        plot(t./60, sigmaVals, ':')
        plot(t./60, -sigmaVals, ':')
        hold off
        ylim(yBounds(sigmaVals, vals));
    end
end


end