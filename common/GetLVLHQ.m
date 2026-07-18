function q = GetLVLHQ(r, v)
%GetLVLHQ Generates the LVLH frame quaternion (ECI -> LVLH)
%   q = GetLVLHQ(r, v)
%
%   Local Vertical Local Horizontal (LVLH) is defined as: z in direction of
%   -r, y in direction of negative orbit normal (-r x v), 
%   x completes triad (y x z)
%   Src: Markley 2014 pg.36
%
%   Inputs
%   r   (3, 1)  Orbital position in ECI
%   v   (3, 1)  Orbital velocity in ECI
%
%   Outputs
%   q   (4, 1)  ECI -> LVLH quaternion

% directions
z = -r;
y = cross(-r, v);

% normalise
y = y./norm(y);
z = z./norm(z);

x = cross(y, z);

% attitude (ECI -> LVLH)
A = [x y z]';
q = DCMToQ(A);

end