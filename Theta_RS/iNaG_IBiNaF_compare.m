function iNaG_IBiNaF_compare

v = -100:50;

[iNaG_m, iNaG_h, iNaG_ht] = iNaG_activation(v);

[iNaF_m, iNaF_h, iNaF_ht] = IBiNaF_activation(v);

figure

subplot(311)

plotyy(v', [iNaG_m; iNaF_m]', v', [iNaG_h; iNaF_h]')

subplot(312)

plotyy(v', [(iNaG_m.^3).*iNaG_h; (iNaF_m.^3).*iNaF_h]', v', [iNaG_ht; iNaF_ht]')

subplot(313)

plotyy(v', [(iNaG_m.^3).*iNaG_h - (iNaF_m.^3).*iNaF_h]', v', [iNaG_ht - iNaF_ht]')

