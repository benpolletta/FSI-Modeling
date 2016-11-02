function [minf, hinf, htau] = iNaG_activation(v)

Noffset = 0;

am = -(v+16-Noffset)./(10*(exp(-(v+16-Noffset)/10)-1));
bm = 4*exp(-(v+41-Noffset)/18);

minf = am./(am + bm);

ah = .07*exp(-(v+30-Noffset)/20);
bh = 1./(exp(-(v-Noffset)/10)+1);

hinf = ah./(ah + bh);

temp = 22;

htau = ones(size(v))/5; % /(3^((temp-6.3)/10));