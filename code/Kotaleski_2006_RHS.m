function [V_out] = Kotaleski_2006_RHS(t,V_in,no_cells,I_app,I_on,noise_multiplier)

    gNa = 114.9;  ENa = 45;             %Sodium max conductance and reversal.

    gK3132 = 58.2;  EK = -90;           %Potassium max conductance and reversal.
    
    gK13 = 0.146;
    
    gKA = 33.3; tau_hKA = 14;
    
    gL = 0.25;  ERest = -70;              %Leak max conductance and reversal.

    Vs_in = V_in(1:no_cells);
    m_in = V_in(no_cells + (1:no_cells));
    h_in = V_in(2*no_cells + (1:no_cells));
    n13_in = V_in(3*no_cells + (1:no_cells));
    n3132_in = V_in(4*no_cells + (1:no_cells));
    mKA_in = V_in(5*no_cells + (1:no_cells));
    hKA_in = V_in(6*no_cells + (1:no_cells));
    
    Vs_out = gNa.*(m_in.^3).*h_in.*(ENa - Vs_in) + gK3132.*(n3132_in.^2).*(EK-Vs_in) ...
      + gK13.*(n13_in.^4).*(EK - Vs_in) + gKA.*(mKA_in.^4).*hKA_in.*(EK - Vs_in) ...
      + gL*(ERest-Vs_in) + I_app*(t > I_on) + noise_multiplier*randn; %real(sqrt(2*noise_multiplier*randn));    %Update I-cell voltage.
    m_out = dx_dt(alpha_m(V), beta_m(V), m_in);
    h_out = dx_dt(alpha_h(V), beta_h(V), h_in);                                            %Update h.
    n13_out = dx_dt(alpha_n13(V), beta_n13(V), n13_in);                                            %Update n.
    n3132_out = dx_dt(alpha_n3132(V), beta_n3132(V), n3132_in);                                                    %Update a.
    mKA_out = (m_inf(Vs_in) - mKA_in)./tau_m(Vs_in);                                                    %Update b.
    hKA_out = (h_inf(Vs_in) - hKA_in)./tau_hKA;
    
    V_out = [Vs_out; m_out; h_out; n13_out; n3132_out; mKA_out; hKA_out];
    
end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function x_out = dx_dt(alpha,beta,x)
x_out = alpha*(1 - x) - beta*x;
end

function a_M = alpha_m(V)
num = 3020 - 40*V;
denom = exp((V+75)/13.5) - 1;
a_M = num/denom;
end

function b_M = beta_m(V)
num = 1.2262;
denom = exp(V/42.248);
b_M = num./denom;
end

function a_H = alpha_h(V)
num = .0035;
denom = exp(V/24.186);
a_H = num./denom;
end

function b_H = beta_h(V)
num = -(0.8712 + 0.017*V);
denom = exp((V + 51.25)/(-5.2)) - 1;
b_H = num./denom;
end

function a_N = alpha_n13(V)
num = -(0.616 + 0.014*V);
denom = exp((V+44)/(-2.3)) - 1;
a_N = num./denom;
end

function b_N = beta_n13(V)
num = 0.0043;
denom = exp((V+44)/(34));
b_N = num./denom;
end

function a_N = alpha_n3132(V)
num = -(95 - V);
denom = exp((V-95)/(-11.8)) - 1;
a_N = num./denom;
end

function b_N = beta_n3132(V)
num = 0.025;
denom = exp(V/22.222);
b_N = num./denom;
end

function m_i = m_inf(V)
num = 0.001;
denom = exp((V + 45)/(-13)) + 1;
m_i = num./denom;
end

function t_M = tau_m(V)
t_M = 0.000001*(exp((V + 70)/(-13)) + 1);
end

function h_i = h_inf(V)
num = 0.001;
denom = exp((V + 77)/8) + 1;
h_i = num./denom;
end