function iKDRG_IBiKDR_compare

v = -100:50;

[iKDRG_minf, iKDRG_mtau] = iKDRG_activation(v);

[IBiKDR_minf, IBiKDR_mtau] = IBiKDR_activation(v, 0);

figure

subplot(211)

[~, h1, h2] = plotyy(v', [iKDRG_minf.^4; IBiKDR_minf.^4]', v', [iKDRG_mtau; IBiKDR_mtau]');

title('Potassium Activation Curves', 'FontSize', 16)

legend(h1, {'Gutfreund''s m_{\infty}^4', 'Dave''s m_{\infty}^4'}, 'Location', 'E', 'FontSize', 12)

legend(h2, {'Gutfreund''s \taum', 'Dave''s \taum'}, 'Location', 'W', 'FontSize', 12)

subplot(212)

[~, h1, h2] = plotyy(v', [iKDRG_minf.^4 - IBiKDR_minf.^4]', v', [iKDRG_mtau - IBiKDR_mtau]');

legend(h1, {'Gutfreund''s m_{\infty}^4 - Dave''s m_{\infty}^4'}, 'Location', 'SE', 'FontSize', 12)

legend(h2, {'Gutfreund''s \taum - Dave''s \taum'}, 'Location', 'W', 'FontSize', 12)

xlabel('Voltage (mV)', 'FontSize', 14)


save_as_pdf(gcf, 'iKDRG_IBiKDR_compare')

