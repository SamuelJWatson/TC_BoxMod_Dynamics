function [SH] = SHsat(Temp)
%SHsat Calculate specific humidity at saturation
%   function calculated using Curve Fitting 'exp2'

a = 1.445e-06;
b = 0.2205; 
c = 4.967;
d = 0.05718;

SH = a*exp(b*Temp) + c*exp(d*Temp);
end