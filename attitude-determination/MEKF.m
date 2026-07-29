classdef(Abstract) MEKF < handle
%MEKF Abstract superclass for MEKF subclasses
%   
%   Subclasses:
%       MEKFNominal
%       MEKFGyroCal
%   
%   Methods:
%       Initialise
%       Step
%       CalcError
%       CalcSigma3
%       Plot
properties(Abstract, Constant)
    nStates (1, 1) double {mustBeInteger, mustBePositive}
end
properties(SetAccess=protected)
    xEst (:, 1) double
    qEst (4, 1) double
    P (:, :) double

    wEst (3, 1) double

    STrue (3, 3) double

    Q (:, :) double
    R (:, :) double

    F (:, :) double
    G (:, :) double
    y (:, 1) double
    h (:, 1) double
    H (:, :) double
    K (:, :) double
end

methods(Static)
    %%% compose the calibration matrix S
    function S = SMatrix(s, kU, kL)
        S = [
            s(1) kU(1) kU(2);
            kL(1) s(2) kU(3);
            kL(2) kL(3) s(3)
        ];
    end

    %%% gyro measurement function
    function wGyro = GyroMeasurement(wTrue, STrue, bTrue, sig)
        wGyro = (eye(3) + STrue)*wTrue + bTrue + sig*randn(3, 1);
    end
    
    %%% discrete covariance propagation
    function PNew = PropPDisc(POld, F, G, Q, nStates, dt)
        A = [
            -F, G*Q*G';
            zeros(nStates, nStates), F';
        ].*dt;
    
        B = expm(A);
    
        Phi = B(nStates+1:end, nStates+1:end)';
        QK = Phi*B(1:nStates, nStates+1:end);
    
        PNew = Phi*POld*Phi' + QK;
    end
    
    %%% discrete quaternion propagation
    function qNew = PropQDisc(qOld, w, dt)
        dw = 0.5*norm(w)*dt;
    
        psi = sin(dw)*w/norm(w);
    
        Theta = [
            cos(dw), -psi';
            psi, cos(dw)*eye(3)-Skew(psi);
        ];
    
        qNew = Theta*qOld;
    end

    function Plot(figData)
        %Plot Creates a tiledlayout figure of state errors for subclasses
        %   Plot(figData)
        %   
        %   Inputs
        %   t           (1, :)      List of timesteps (s)
        %   xErrList    (n, :)      List of error states
        %   sigma3List  (n, :)      List of 3-sigma bounds
        %   figData     (struct)    Information to create the figure
        %   .figName    (string)    Name of figure
        %   .tlDim      (2, 1)      Number of [rows, cols] of tiledlayout
        %   .plotSizes  (n/3, 2)    Number of tiles that each plot takes up
        %                           [1 1] for 1 row, 1 col
        %                           [1 2] for 1 row, 2 col etc.
        %   .yLabels    (n/3, 1)    String array of plot y-axis labels
        %   .xLabel     string      x-axis label for all plots
        %   .y          {n/3, 1}    Cell array of state errors to plot
        %   .sigma3s    {n/3, 1}    Cell array of 3-sigma bounds to plot
        %   .x          (:, 1)      Timesteps to plot

        figure("Name", figData.figName);

        tl = tiledlayout(figData.tlDim(1), figData.tlDim(2));
        tl.Title.String = "MEKF State Errors";
        colororder(lines(3))
        set(0, 'DefaultLineLineWidth', 1.4)

        nPlots = size(figData.plotSizes, 1);

        for i = 1:nPlots
            nexttile(figData.plotSizes(i, :))
            plot(figData.x, figData.y{i})
            grid on
            ylabel(figData.yLabels(i))
            xlabel(figData.xLabel)
            hold on
            plot(figData.x, figData.sigma3s{i}, ':')
            plot(figData.x, -figData.sigma3s{i}, ':')
            hold off
            ylim(yBounds(figData.sigma3s{i}, figData.y{i}));

            if (i == 1)
                lg = legend('x', 'y', 'z', '3-\sigma_x', '3-\sigma_y', '3-\sigma_z');
                lg.Location = "eastoutside";
            end
        end

        function bound = yBounds(sigVals, xVals)
            % sets smart y limits
            
             % get the timestep with smallest sigma-3 bound
            [~, minI] = min(sum(sigVals));

            % get max of [sigma-3; x] at this timestep
            % set y limits to 3 times this value
            y = 3*(max([sigVals(:, minI); xVals(:, minI)]));

            bound = [-y y];
        end
    end
end
methods(Abstract, Access=protected)
    wEst = wEstFn(o, wGyro)
    [y, h, H] = MeasurementMatrices(o, extData)
    [F, G] = StateMatrices(o)
end
methods
    function o = Initialise(o, x0, q0, P0, STrue, Q, R)
        %Initialise Initialises the filter with initial values
        %   o = Initialise(o, x0, q0, P0, STrue, Q, R)
        %   
        %   Inputs
        %   x0      (n, 1)  Initial estimated state vector
        %   q0      (4, 1)  Initial estimated quaternion
        %   P0      (n, n)  Initial covariance matrix
        %   STrue   (3, 3)  True gyroscope calibration matrix
        %   Q       (n, n)  Spectral density matrix
        %   R       (:, :)  Measurement covariance matrix

        o.xEst = x0;
        o.qEst = q0;
        o.P = P0;
        o.STrue = STrue;
        o.Q = Q;
        o.R = R;
    end
    
    function o = Step(o, dt, extData)
        %Step Steps the filter forward by one timestep
        %   o = Step(o, dt, extData)
        %   
        %   Use this function within the simulation loop
        %
        %   Inputs
        %   dt      (1, 1)      Loop timestep (s)
        %   extData (struct)    Data struct of values at this timestep
        %   .wGyro  (3, 1)      Measured gyro angular velocity (rad/s)
        %   .~                  Other fields as required by subclasses
        %
        %   Output
        %   o       (object)    MEKF object instance

        % simulate measurement
        [o.y, o.h, o.H] = o.MeasurementMatrices(extData);
    
        % calculate gain
        o.K = o.P*o.H' / (o.H*o.P*o.H' + o.R);
    
        % measurement update
        % covariance update
        o.P = (eye(o.nStates)-o.K*o.H)*o.P;
        % state update
        o.xEst(1:3) = zeros(3, 1);
        o.xEst = o.xEst + o.K*(o.y - o.h);
        % quaternion update
        o.qEst = o.qEst + QKinematics([o.qEst; o.xEst(1:3)], struct('iQ', 1:4, 'iWSat', 5:7));
        o.qEst = o.qEst./norm(o.qEst);
    
        % propagate dynamics
        o.wEst = o.wEstFn(extData.wGyro);
        [o.F, o.G] = o.StateMatrices();
    
        o.P = o.PropPDisc(o.P, o.F, o.G, o.Q, o.nStates, dt);
        o.qEst = o.PropQDisc(o.qEst, o.wEst, dt);
    end

    function xErr = CalcError(o, qTrue, xiTrue)
        %CalcError Generates the state error
        %   xErr = CalcError(o, qTrue, xiTrue)
        %
        %   Inputs
        %   qTrue   (4, 1)      True quaternion
        %   xiTrue  (n-3, 1)    True states EXCLUDING attitude state
        %
        %   Outputs
        %   xErr    (n, 1)      Error state vector
        
        xErr = zeros(o.nStates, 1);

        qErr = QProd(qTrue, QConj(o.qEst));
        qErr = qErr/norm(qErr); 
        angErr = 2*qErr(2:4)/qErr(1);

        xErr(1:3) = angErr;
        xErr(4:end) = xiTrue - o.xEst(4:end);
    end

    function xSigma3 = CalcSigma3(o)
        %CalcSigma3 Calculates the 3-sigma bounds
        %   xSigma3 = CalcSigma3(o)
        %
        %   Outputs
        %   xSigma3    (n, 1)      3-sigma bounds for each state
        
        xSigma3 = 3*sqrt(max(diag(o.P),0));
    end
end
end