function [Ang, Amp, Cmplx] = func_remapAng (in)
% func_remapAng recievies a complex map, ans shifts its phase and returns phase map without changing its amplitude
% this function is used because matlab angle operation always returns values from -pi to pi

Ang = angle(in);
angind = Ang < -pi/2;
Ang(angind) = Ang(angind) + 2*pi;
Ang = Ang;% - pi;

Amp = abs(in);

RealX = Amp.*cos(Ang);
ImaqY = Amp.*sin(Ang);
clear i
Cmplx = RealX + 1i * ImaqY;