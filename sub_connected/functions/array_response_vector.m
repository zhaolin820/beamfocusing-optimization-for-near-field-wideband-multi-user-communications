function [a] = array_response_vector(r, theta, N, d, f)
%The near-field array response vector for ULAs
%  [a] = array_response_vector(r, theta, N, d, f)
%Inputs:
%   r: distance
%   theta: angle of direction
%   N: number of antennas
%   d: antenna spacing
%   f: carrier frequency
%Outputs:
%   a: array response vector
%Date: 04/04/2024
%Author: Zhaolin Wang

c = 3e8;

n = 0:(N-1);
delta_n = (n-(N-1)/2) * d;
delta_n = delta_n';

% distance
r_n = sqrt(r^2 + delta_n.^2 - 2*r*delta_n*cos(theta));

% array response vector
a = exp( -1i * 2 * pi * f/c * r_n );

end

