classdef MEKFNominal < MEKF
%MEKFNominal A 6-state MEKF involving attitude and gyro biases
%   mekf = MEKFNominal
%
%   State Vector
%   x (6) = [
%           attitude error (3); 
%           gyro biases (3)
%   ]
%
%   See the MEKF parent class for more documentation on methods

properties(Constant)
    nStates = 6;
end
properties 
    C (:, :) double
end

methods
    % add the matrix C to initialisation
    function o = Initialise(o, x0, q0, P0, STrue, Q, R)
        Initialise@MEKF(o, x0, q0, P0, STrue, Q, R);

        o.C = eye(3) / (eye(3) + STrue);
    end
end
methods(Access=protected)
    function wEst = wEstFn(o, wGyro)
        bEst = o.xEst(4:6);
        
        wEst = o.C*(wGyro - bEst);
    end
    % TODO: once magnetometer model is written, use y directly from
    % that
    function [y, h, H] = MeasurementMatrices(o, extData)
        dMag = MeasMagnetometerEarth;
        dMag.jD = extData.jD;
        dMag.kR = 1:3;
        dMag.kQ = 4:7;
        dMag.quantization = 1e-20;
        x = [extData.r; extData.qTrue;];
    
        bMeas = MagModel(x, dMag);
        % TODO: use our own magnetic reference function (should be done
        % when mag model is done)
        bEst = QToDCM(o.qEst)*BDipole(extData.r, extData.jD);
    
        % TODO: add sigmaMag (and sigmaSun eventually) to class
        % properties
        y = bMeas + extData.sigmaMag*randn(3, 1);
    

        H = [Skew(bEst) zeros(3, o.nStates-3)];
        h = bEst;
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
        %   t           (1, :)  List of timesteps (s)
        %   xErrList    (n, :)  List of error states
        %   sigma3List  (n, :)  List of 3-sigma bounds
        
        figure;
        tl = tiledlayout(2, 1);
        tl.Title.String = "MEKF State Errors";
        colororder(lines(3))
        set(0, 'DefaultLineLineWidth', 1.4)
        
        nexttile
        plot(t./60, rad2deg(xErrList(1:3,:)))
        grid on
        ylabel('\delta\vartheta (deg)')
        xlabel('t (min)')
        hold on
        plot(t./60, rad2deg(sigma3List(1:3, :)), ':')
        plot(t./60, -rad2deg(sigma3List(1:3, :)), ':')
        hold off
        
        lg = legend('x', 'y', 'z', '3-\sigma_x', '3-\sigma_y', '3-\sigma_z');
        lg.Location = "eastoutside";
        
        nexttile
        plot(t./60, rad2deg(xErrList(4:6, :)).*3600)
        grid on
        ylabel('\Delta\beta (deg/h)')
        xlabel('t (min)')
        hold on
        plot(t./60, rad2deg(sigma3List(4:6, :)).*3600, ':')
        plot(t./60, -rad2deg(sigma3List(4:6, :)).*3600, ':')
        hold off
    end
end
end