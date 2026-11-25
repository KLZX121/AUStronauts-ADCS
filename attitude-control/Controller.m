classdef Controller < handle

properties

    % minimum angular velocity (deg/s) for detumbling to meet
    detumbleThreshold double
    
    % minimum angular error (deg) for pointing to meet
    pointThreshold double

    % mission phases
    missionSequence {}
    nPhases int16 = 0
    % current index of mission phase
    currentPhase int16 = 1
    % timestep that current phase was completed
    t0 double = -1
    % whether the phase has been completed (to start counting t0)
    phaseCompleted logical = false

    % whether to print logging information
    doLogging logical = false

    psiPrev double = 0
end
methods
    function o = Controller(missionSequence, detumbleThreshold, pointThreshold, doLogging)
        %   o = Controller(missionSequence, detumbleThreshold, pointThreshold, doLogging)
        %   
        %   Initiates the controller with the mission phases
        %
        %   The input is given as a {:, 3} cell array in the form:
        %   {"mode1", "target", t1;
        %    "mode2", "target", t2; ...}
        %   where 
        %   "mode1" and "mode2" are consecutive control modes defined
        %       in the <currentMode> property
        %       ["none","detumble","point"]
        %   "target" is what pointing mode is desired (n/a for detumble)
        %       ["nadir","sun","hybrid"]
        %   +t1 and +t2 are the amount of timesteps to hold the control
        %       mode for counted AFTER the threshold accuracy is achieved.
        %       note that it is simulation step, not actual time
        %   
        %   Inputs
        %   missionSequence     {:, 3}  Cell array
        %   detumbleThreshold   (1, 1)  Maximum angular velocity (deg) for
        %                               detumbling to count as completed
        %   pointThreshold      (1, 1)  Max angular error (deg) for
        %                               pointing to count as completed
        %   doLogging           bool    Whether to print logs to console

        o.missionSequence = missionSequence;
        o.nPhases = size(missionSequence, 1);

        o.detumbleThreshold = detumbleThreshold;
        o.pointThreshold = pointThreshold;

        o.doLogging = doLogging;
    end
    function [TRW, TMTQ, info] = Run(o, step, x, d, env)
    %Run    Runs the controller at each timestep   
    %   [TRW, TMTQ, info] = Run(o, step, x, d, env)
    %
    %   Inputs
    %   step    (1, 1)  Current simulation step
    %   x       (:, 1)  State Vector
    %   d       struct  Data struct
    %   env     struct  Environmental struct
    %   .B      (3, 1)  Local magnetic field (ECI) (T)
    %   .s      (3, 1)  Sun direction vector (sat -> sun) ECI
    %
    %   Outputs
    %   TRW     (3, 1)  Reaction wheel torque (internal) (N m)
    %   TMTQ    (3, 1)  Magnetorquer torque (external) (N m)
    %   info    struct  Supplementary info
    %   .m      (3, 1)  Magnetorquer dipole moment (A m^2)

        if o.t0 == -1
            o.t0 = step;
            if o.doLogging
                disp('Controller initialised')
                fprintf('Current mode: %s %s\n', o.missionSequence{o.currentPhase, 2}, o.missionSequence{o.currentPhase, 1})
            end
        end

        if o.currentPhase > o.nPhases
            o.currentPhase = o.nPhases;
        end
    
        % get current phase info
        [phaseMode, phaseTarget, phaseT] = o.missionSequence{o.currentPhase, :};
    
        TRW = zeros(3, 1);
        TMTQ = zeros(3, 1);
        info.m = zeros(3, 1);
        if phaseMode == "none"
            % apply no torques
            if ~o.phaseCompleted
                o.t0 = step;
                o.phaseCompleted = true;
            end
        elseif phaseMode == "detumble"
            % calculate detumbling torque
            [TMTQ, info.m] = o.DetumbleTorque(x, d, env);
    
            currentAngVelDeg = rad2deg(norm(x(d.iWSat)));

            if ~o.phaseCompleted && (currentAngVelDeg <= o.detumbleThreshold)
                o.t0 = step;
                o.phaseCompleted = true;

                if o.doLogging
                    fprintf('Detumbling completed at step %d\n', step)
                    fprintf('Angular velocity %.4f deg/s vs threshold of %.4f deg/s\n', currentAngVelDeg, o.detumbleThreshold)
                end
            end
        elseif phaseMode == "point"
            % calculate desired attitude based on target
            qTarget = o.CalcDesiredQuat(phaseTarget, x, d, env);

            % calculate control torque
            [TC, degErr] = o.ControlTorque(qTarget, x, d);
            TRW = -TC;
    
            if ~o.phaseCompleted && (degErr <= o.pointThreshold)
                o.t0 = step;
                o.phaseCompleted = true;

                if o.doLogging
                    fprintf('%s pointing completed at step %d\n', phaseTarget, step)
                    fprintf('Angular error %.4f deg vs threshold of %.4f deg\n', degErr, o.pointThreshold)
                end
            end
        end
    
        % if the +t amount of time has passed, move onto the next phase
        if o.phaseCompleted && ((step - o.t0) >= phaseT)
           o.currentPhase = min(o.currentPhase + 1, o.nPhases);
           o.phaseCompleted = false;

           if o.doLogging
                fprintf('Moving to next mode: %s %s at step %d\n', o.missionSequence{o.currentPhase, 2}, o.missionSequence{o.currentPhase, 1}, step)
           end
        end
    end

% ======================== CONTROL LAWS ========================

    function [T, m] = DetumbleTorque(~, x, d, env)
        %   INPUTS
        %   x               (:, 1)  State vector
        %   d               struct  Input parameters
        %   .iWSat          (1, 3)  Indices of satellite ang vel in state vector
        %   env             struct  Environmental Data
        %   .B              (3, 1)  Local magnetic field (ECI) (T)
        %
        %   OUTPUT
        %   T               (3, 1)  Control torque vector (x,y,z) (N m)
        %   m               (3, 1)  Magnetorquer dipole moment (A m^2)

        w = x(d.iWSat); % current satellite angular velocities

        % Gain for controller
        k = 4*pi/(90*60)*(1+sin(51.6*pi/180))*0.0067;
        m = cross(k/norm(env.B)*w,env.B/norm(env.B));
        T = cross(m,env.B);
    end
    function [T, deg_q] = ControlTorque(~, q_d, x, d)
        %   INPUTS
        %   q_d             (4, 1)  desired direction quaternion input
        %   x               (:, 1)  state vector
        %   d                       Input parameters (struct)
        %       .iWSat      (1, 3)  Indices of satellite ang vel in state vector
        %       .iWRW       (1, 3)  Indices of rw ang vel in state vector
        %
        %   OUTPUT
        %   T               (3, 1)  control torque vector (x,y,z)
        %   deg_q           (1, 1)  angular error between current and
        %                           desired attitude (deg)

        % Find delta quaternion

        % input
        q_current = x(d.iQ);
        q_desired = q_d;

        % Find the error quaternion (change in quaternion needed)
        delta_q = QProd( q_current, QConj(q_desired));
        w = x(d.iWSat); % current satellite angular velocities

        %temporary fix, change input to include a mode (slew or track)
        % code to find amount of turn
        deg_q = 2*acos(abs(delta_q(1)/norm(delta_q(1:4))))*180/pi;
        if deg_q < 5
            operating_mode = "small_delq";
        else
            % Big delta q (large angle flip, move to LVLH, slew)
            operating_mode = "slew";
        end

        if operating_mode == "small_delq"

            k_p = 1.62e-4; % proportional gain (CHANGE)
            k_d = 9.75e-4; % derivative gain (CHANGE)
            % Calculate the control torque vector
            T = -k_p * sign(delta_q(1)) * delta_q(2:4) - k_d ...
                * (1 + delta_q(2:4).' * delta_q(2:4)) * w;

        elseif operating_mode == "slew"

            % Calculate the control torque vector

            %step 1: find angle between quaternions: deg_q
            %step 2: choose interpolation coefficient: T
            T = 0.82/deg_q;

            %step 3: interpolate to interpolation coefficient using slerp
            q_int = sin((1-T)*deg_q*pi/180) / sin(deg_q*pi/180) * -q_current ...
                + sin(T*deg_q*pi/180) / sin(deg_q*pi/180) * q_desired;

            %step 4: find error quaternion to interpolated quaternion
            delta_q_int = QProd( q_current, QConj(q_int));

            %step 5: calculate the control torque vector
            k_p = 1.62e-4; % proportional gain (CHANGE)
            k_d = 9.75e-4; % derivative gain (CHANGE)
            T = -k_p * delta_q_int(2:4) - k_d * w;

            %deg_q_int = 2*acos(abs(delta_q_int(1)/norm(delta_q_int(1:4))))*180/pi %debug
        end

        % Limit to max RW speed
        for j = 1:3
            if abs(x(d.iWRW(j))) > d.maxWRW
                T(j) = 0;
            end
        end
    end

% ======================== DESIRED TARGETS ========================

    function q = CalcDesiredQuat(~, target, x, d, env)
    %   Calculates the desired attitude given the target
        if target == "nadir"
            q = GetLVLHQ(x(d.iR), x(d.iV));
    
        elseif target == "sun"
            % aligns the satellite's x axis with the sun direction
            dcm = AlignVecs(env.s, [1; 0; 0]);
            q = DCMToQ(dcm);
        
        elseif target == "hybrid"
            % aligns roll and pitch with lvlh but yaw with sun
            qN = GetLVLHQ(x(d.iR), x(d.iV));

            % sun direction in LVLH axes
            s = env.s(:) / norm(env.s);
            sL = QToDCM(qN)*s;

            % yaw angle of the sun projection measured from LVLH x.
            % If the sun is close to the nadir/zenith axis the projection
            % is ill-conditioned, so hold the previous yaw instead
            if norm(sL(1:2)) < sin(deg2rad(5))
                psi = o.psiPrev;
            else
                psi = atan2(sL(2), sL(1));
            end
            o.psiPrev = psi;

            % frame rotation about z by psi, then compose: A(q) = A(qYaw)*A(qN)
            qYaw = [cos(psi/2); 0; 0; sin(psi/2)];
            q = QProd(qYaw, qN);
            q = q(:) / norm(q);
        end
    end
end
end