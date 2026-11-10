classdef RIEKF < MEKF
%RIEKF A 6-state right-invariant EKF (attitude + gyro bias)
%   riekf = RIEKF
%
%   Based on: Gui & de Ruiter, "Quaternion Invariant Extended Kalman
%   Filtering for Spacecraft Attitude Estimation", JGCD 41(4), 2018.
%
%   State Vector
%   xEst (6) = [
%           attitude error (3);     (always zero, attitude lives in qEst)
%           bias estimate (3)       (body frame, same as MEKF)
%   ]
%
%   The covariance P (and therefore CalcSigma3) is for the INERTIAL-frame
%   error vector [gammaI; betaI], where (paper eq. 14, Hamilton product)
%       exp_q(gammaI/2) = qEst (x) qTrue*
%       betaI           = A'(qTrue)*(bEst - bTrue)
%   For isotropic initial covariances P0 is identical to the MEKF's
%   (paper eq. 40).
%
%   Quaternion convention
%   The attitude quaternion is the same quantity as in the MEKF
%   (scalar first, A(q) maps inertial -> body). Only the quaternion
%   product order differs: the paper uses the Hamilton product
%       q (x) p = [q0*p0 - qv'*pv; q0*pv + p0*qv + cross(qv, pv)]
%   All quaternion maths in this class uses RIEKF.HamProd so it does not
%   depend on the ordering used by QProd.
%
%   Numerical propagation
%   Uses the same exact discretisation as MEKF (closed-form quaternion
%   exponential and Van Loan matrix exponential). This is exact here
%   because, with zero-order-hold gyro, the inertial angular velocity
%   estimate A'(q)*(wGyro - b) is constant over a step: the attitude
%   rotates about w, which leaves w unchanged. F is therefore constant.
%
%   See MEKF superclass for method documentation

properties(Constant)
    nStates = 6;
end

methods(Access=protected)
    function wEst = wEstFn(o, wGyro)
        wEst = wGyro - o.xEst(4:6);
    end

    function [F, G] = StateMatrices(o)
        % paper eq. 27
        A = RIEKF.DCM(o.qEst);
        wI = A'*o.wEst;     % inertial angular velocity estimate

        F = [
            zeros(3), -eye(3);
            zeros(3), Skew(wI);
        ];
        G = blkdiag(A', -A');
    end
end
methods
    function o = Update(o, q, env)
        %Update Right-invariant measurement update (paper eq. 20, 33)
        %   o = Update(o, q, env)
        %
        %   See MEKF.Update for inputs

        % simulate measurement (H is in the inertial-error coordinates)
        [o.y, o.h, o.H, o.R] = o.MeasurementMatrices(q, env);

        % calculate gain
        o.K = o.P*o.H' / (o.H*o.P*o.H' + o.R);

        % covariance update
        o.P = (eye(o.nStates) - o.K*o.H)*o.P;

        % correction in the inertial frame
        c = o.K*(o.y - o.h);
        cq = c(1:3);
        cb = c(4:6);

        % quaternion update: q+ = exp_q(-cq/2) (x) q-
        o.qEst = RIEKF.HamProd(RIEKF.QExp(-cq/2), o.qEst);
        o.qEst = o.qEst./norm(o.qEst);

        % bias update: b+ = b- - A(q+)*cb
        o.xEst(1:3) = zeros(3, 1);
        o.xEst(4:6) = o.xEst(4:6) - RIEKF.DCM(o.qEst)*cb;
    end

    function [y, h, H, R] = MeasurementMatrices(o, q, env)
        %MeasurementMatrices Measurements with H for the inertial error
        %   [y, h, H, R] = MeasurementMatrices(o, q, env)
        %
        %   y, h, R are identical to MEKF.MeasurementMatrices (body-frame
        %   residual y - h). H is re-expressed so that
        %       y - h = H*[gammaI; betaI]
        %   using the MEKF relation y - h = H_B*[dtheta_B; dbeta] with
        %       dtheta_B = -A*gammaI,  dbeta = -A*betaI
        %   i.e. H = -H_B*blkdiag(A, A).
        %
        %   For vector sensors this is equivalent to the paper's inertial
        %   output error E = r - A'(q)*y with H = [r x, 0] (eq. 29-31):
        %   rotating the residual by -A' leaves K*(residual) unchanged.
        %   Working in the body frame lets the scalar CSS measurements be
        %   used as well.

        [y, h, H, R] = MeasurementMatrices@MEKF(o, q, env);

        A = RIEKF.DCM(o.qEst);
        H = -H*blkdiag(A, A);
    end

    function xErr = CalcError(o, qTrue, xiTrue)
        %CalcError Generates the inertial-frame state error
        %   xErr = CalcError(o, qTrue, xiTrue)
        %
        %   Inputs
        %   qTrue   (4, 1)      True quaternion
        %   xiTrue  (3, 1)      True gyro bias (body frame)
        %
        %   Outputs
        %   xErr    (6, 1)      [gammaI; betaI] (estimate - truth, paper
        %                       eq. 14)

        xErr = zeros(o.nStates, 1);

        qErr = RIEKF.HamProd(o.qEst, QConj(qTrue));
        qErr = qErr/norm(qErr);

        xErr(1:3) = RIEKF.QLog(qErr);
        xErr(4:6) = RIEKF.DCM(qTrue)'*(o.xEst(4:6) - xiTrue(1:3));
    end
end
methods(Static)
    function A = DCM(q)
        %DCM Rotation matrix A(q), inertial -> body (paper eq. 1)
        qv = q(2:4);
        A = (q(1)^2 - qv'*qv)*eye(3) + 2*(qv*qv') - 2*q(1)*Skew(qv);
    end

    function r = HamProd(q, p)
        %HamProd Hamilton quaternion product q (x) p (paper convention)
        qv = q(2:4);
        pv = p(2:4);
        r = [
            q(1)*p(1) - qv'*pv;
            q(1)*pv + p(1)*qv + cross(qv, pv);
        ];
    end

    function q = QExp(x)
        %QExp Quaternion exponential map R^3 -> S^3 (paper eq. 4)
        n = norm(x);
        if n < 1e-8
            s = 1 - n^2/6;
        else
            s = sin(n)/n;
        end
        q = [cos(n); s*x];
    end

    function g = QLog(q)
        %QLog Principal rotation vector g such that exp_q(g/2) = q
        if q(1) < 0
            q = -q;
        end
        n = norm(q(2:4));
        if n < 1e-12
            g = 2*q(2:4);
        else
            g = 2*atan2(n, q(1))*q(2:4)/n;
        end
    end

    function Plot(t, xErrList, sigma3List, plotType)
        %Plot Plots the state errors over time, as well as 3-sigma bounds
        %   Plot(t, xErrList, sigma3List, plotType)
        %
        %   Inputs
        %   t           (1, :)      List of timesteps (min)
        %   xErrList    (6, :)      List of error states
        %   sigma3List  (6, :)      List of 3-sigma bounds
        %   plotType    (str)       'overview' or 'detailed'

        yLabels = [
            "\delta\vartheta_I (deg)";
            "\Delta\beta_I (deg/h)"
        ];
        y = {
            rad2deg(xErrList(1:3, :));
            rad2deg(xErrList(4:6, :)).*3600;
        };
        sig3 = {
            rad2deg(sigma3List(1:3, :));
            rad2deg(sigma3List(4:6, :)).*3600;
        };

        if plotType == "overview"
            figData.figName = "RIEKF";
            figData.plotType = plotType;

            figData.tlDim = [2, 1];
            figData.plotSizes = repmat([1 1], 2, 1);
            figData.yLabels = yLabels;
            figData.y = y;
            figData.sigma3s = sig3;

            figData.xLabel = "t (min)";
            figData.x = t;

            MEKF.Plot(figData)
        elseif plotType == "detailed"
            figNames = [
                "RIEKF Attitude Error";
                "RIEKF Gyro Bias Error";
            ];

            for i = 1:length(figNames)
                figData.figName = figNames(i);
                figData.plotType = plotType;

                figData.tlDim = [3, 1];
                figData.plotSizes = repmat([1 1], 3, 1);
                figData.yLabels = [
                    sprintf("X %s", yLabels(i));
                    sprintf("Y %s", yLabels(i));
                    sprintf("Z %s", yLabels(i));
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