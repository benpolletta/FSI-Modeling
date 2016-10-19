function [m_inf, tau_m] = iKs_activation_curves(X)

halfKs = -35; temp = 34;
tau_div = 3^((temp - 22)/10);
tau_mult = 1;

m_inf = 1./(1+exp(-(X - halfKs)/10));
tau_m = tau_mult*1000./(3.3*tau_div*(exp((X - halfKs)/40) + exp(-(X - halfKs)/20)));