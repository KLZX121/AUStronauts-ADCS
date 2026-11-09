function [theta_err, q_err] = QAttErr(q_true, q_est)
%QAttErr Quaternion Attitude Error
%   [theta_err, q_err] = QAttErr(q_true, q_est)
%
%   Computes the attitude error between a true and estimated attitude in
%   quaternion form
%
%   Inputs
%   q_true  (4, 1)  True attitude quaternion
%   q_est   (4, 1)  Estimated attitude quaternion
%
%   Outputs
%   theta_err       Pointing error (rad)
%   q_err   (4, 1)  Error quaternion

q_err = QProd(q_true, QConj(q_est));

% keep sign consistent to prevent flips
if (q_err(1) < 0)
    q_err = q_err.*-1;
end

% renormalise
q_err = q_err/norm(q_err);

% pointing error
theta_err = 2*acos(q_err(1));

end