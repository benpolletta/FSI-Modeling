function [minf, mtau] = IBiKDR_activation(X, KDR_offset)

KDR_V1 = [29.5];
KDR_d1 = [10];
KDR_V2 = [10];
KDR_d2 = [10];

minf = 1./(1+exp((-X-KDR_V1+KDR_offset)/KDR_d1));
mtau = .25+4.35*exp(-abs(X+KDR_V2-KDR_offset)/KDR_d2);
 
