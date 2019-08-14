function r = fmcw_phase2range(phi,lambdac,rc,K,ci)

% r = fmcw_phase2range(phi,lambdac,rc,K,ci)
%
% Convert phase difference to range for FMCW radar
%
% args:
% phi: phase (radians), must be of spectrum after bin centre correction
% lambdac: wavelength (m) at centre frequency
% 
% optional args: (used for precise method)
% rc: coarse range of bin centre (m)
% K = chirp gradient (rad/s/s)
% ci = propagation velocity (m/s)
%
% Craig Stewart
% 2014/6/10

if nargin<5
    % First order method
    r = lambdac*phi./(4*pi);
else
    % Precise
    r = phi./((4*pi/lambdac) - (4*rc*K/ci^2));
end