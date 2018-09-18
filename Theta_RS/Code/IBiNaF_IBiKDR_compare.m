function IBiNaF_IBiKDR_compare(N_offset, K_offset)

v = -100:50;

[iNaF_m, iNaF_h, iNaF_ht] = IBiNaF_activation(v, N_offset);

[IBiKDR_minf, IBiKDR_mtau] = IBiKDR_activation(v, K_offset);

figure

subplot(211)

[ax, ~, ~] = plotyy(v', ((iNaF_m.^3).*iNaF_h)', v', (IBiKDR_minf.^4)');

ylabel(ax(1), 'Na m^3h')

ylabel(ax(2),'K m_{\infty}^4')

title('Activation Curves')

subplot(313)

plot(v', [iNaF_ht; IBiKDR_mtau]');

legend({'Na \tau_h', 'K \tau_m'}, 'Location', 'NW', 'FontSize', 12)

title('Time Constants')

xlabel('Voltage (mV)', 'FontSize', 14)


save_as_pdf(gcf, sprintf('IBiNaF_IBiKDR_compare_Noff%.3g_Koff%.3g', N_offset, K_offset))

