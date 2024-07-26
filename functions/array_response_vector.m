function [a] = array_response_vector(r, theta, N, d, f)
%To calculate the near-field array response of an ULA at a given location
%  [a] = array_response_vector(r, theta, N, d, f)
%Inputs:
%   r: distance
%   theta: angle
%   N: number of antennas of the ULA
%   d: antenna spacing of the ULA
%   f: carrier frequency
%Outputs:
%   a: near-field array response vector
%Date: 22/07/2024
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

