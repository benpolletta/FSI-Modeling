function iKDRG_IBiKDR_compare

v = -100:50;

[iKDRG_minf, iKDRG_mtau] = iKDRG_activation(v);

[IBiKDR_minf, IBiKDR_mtau] = IBiKDR_activation(v);

figure

subplot(211)

plotyy(v', [iKDRG_minf.^4; IBiKDR_minf.^4]', v', [iKDRG_mtau; IBiKDR_mtau]')

subplot(212)

plotyy(v', [iKDRG_minf.^4 - IBiKDR_minf.^4]', v', [iKDRG_mtau - IBiKDR_mtau]')

