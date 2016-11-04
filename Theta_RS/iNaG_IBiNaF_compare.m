function iNaG_IBiNaF_compare

v = -100:50;

[iNaG_m, iNaG_h, iNaG_ht] = iNaG_activation(v);

[iNaF_m, iNaF_h, iNaF_ht] = IBiNaF_activation(v);

figure

subplot(311)

[~, h1, h2] = plotyy(v', [iNaG_m; iNaF_m]', v', [iNaG_h; iNaF_h]');

title('Sodium Activation Curves', 'FontSize', 16)

legend(h1, {'Gutfreund''s m_{\infty}', 'Dave''s m_{\infty}'}, 'Location', 'E', 'FontSize', 12)

legend(h2, {'Gutfreund''s h_{\infty}', 'Dave''s h_{\infty}'}, 'Location', 'W', 'FontSize', 12)

subplot(312)

[~, h1, h2] = plotyy(v', [(iNaG_m.^3).*iNaG_h; (iNaF_m.^3).*iNaF_h]', v', [iNaG_ht; iNaF_ht]');

legend(h1, {'Gutfreund''s m^3h', 'Dave''s m^3h'}, 'Location', 'E', 'FontSize', 12)

legend(h2, {'Gutfreund''s \tau_h', 'Dave''s \tau_h'}, 'Location', 'W', 'FontSize', 12)

subplot(313)

[~, h1, h2] = plotyy(v', [(iNaG_m.^3).*iNaG_h - (iNaF_m.^3).*iNaF_h]', v', [iNaG_ht - iNaF_ht]');

legend(h1, {'Gutfreund''s m^3h - Dave''s m^3h'}, 'Location', 'NW', 'FontSize', 12)

legend(h2, {'Gutfreund''s \tau_h - Dave''s \tau_h'}, 'Location', 'SE', 'FontSize', 12)

xlabel('Voltage (mV)', 'FontSize', 14)


save_as_pdf(gcf, 'iNaG_IBiNaF_compare')

