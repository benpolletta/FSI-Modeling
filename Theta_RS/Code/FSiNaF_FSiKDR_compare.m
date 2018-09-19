function FSiNaF_FSiKDR_compare(N_offset, K_offset)

v = -100:50;

[iNaF_m, iNaF_h, iNaF_ht] = FSiNaF_activation(v, N_offset);

[FSiKDR_minf, FSiKDR_mtau] = FSiKDR_activation(v, K_offset);

figure

subplot(211)

[ax, ~, ~] = plotyy(v', ((iNaF_m.^3).*iNaF_h)', v', (FSiKDR_minf.^4)');

ylabel(ax(1), 'Na m^3h')

ylabel(ax(2),'K m_{\infty}^4')

title('Activation Curves')

subplot(313)

plot(v', [iNaF_ht; FSiKDR_mtau]');

legend({'Na \tauh', 'K \taum'}, 'Location', 'NW', 'FontSize', 12)

title('Time Constants')

xlabel('Voltage (mV)', 'FontSize', 14)


save_as_pdf(gcf, sprintf('FSiNaF_FSiKDR_compare_Noff%.3g_Koff%.3g', N_offset, K_offset))

