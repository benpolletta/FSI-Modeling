function iKs_IBiMMich_compare

X = (-100:50)';

[MM_inf, tau_MM] = IBiMMich_activation_curves(X);

[iKs_inf, tau_iKs] = iKs_activation_curves(X);

figure

subplot(211)

plotyy(X, [iKs_inf MM_inf], X, [tau_iKs tau_MM])

subplot(212)

plotyy(X, iKs_inf - MM_inf, X, tau_iKs - tau_MM)