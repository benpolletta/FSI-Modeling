function [m_inf, taum] = iKs_activation_curves(X)

halfKs = -35; temp = 34;
tau_div = 3^((temp - 22)/10);
taumult = 1;

m_inf = 1./(1+exp(-(X - halfKs)/10));
taum = taumult*1000./(3.3*tau_div*(exp((X - halfKs)/40) + exp(-(X - halfKs)/20)));