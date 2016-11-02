function [minf, hinf, htau] = IBiNaF_activation(v)

NaF_V0 = [34.5];
NaF_V1 = [59.4];
NaF_d1 = [10.7];
NaF_V2 = [33.5];
NaF_d2 = [15];
NaF_c0 = [0.15];
NaF_c1 = [1.15];
NaF_offset = 15;

hinf = 1./(1+exp((v+NaF_V1-NaF_offset)/NaF_d1));

htau = NaF_c0 + NaF_c1./(1+exp((v+NaF_V2-NaF_offset)/NaF_d2));

minf = 1./(1+exp((-v-NaF_V0+NaF_offset)/10));