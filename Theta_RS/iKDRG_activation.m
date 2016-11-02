function [minf, mtau] = iKDRG_activation(X)

temp = 22;
mtau = ones(size(X))/5; % (3^((temp-6.3)/10));
Koffset = 0;

% Functions
aM = -.01*(X+20-Koffset)./(exp(-(X+20-Koffset)/10)-1);
bM = .125*exp(-(X+30-Koffset)/80);

minf = aM./(aM + bM);
