function [V_out] = Golomb_2007_gating_vars(t,V_in,CS,CG,g_struct,I_app,I_on,noise_multiplier)

    ENa=50;                             %Sodium reversal.
    sigma_m = 11.5;
    theta_h = -58.3; sigma_h = -6.7;

    EK=-90;                             %Potassium reversal.
    theta_n = -12.4; sigma_n = 6.8;     %Delayed rectifier parameters.

    theta_a = -50; sigma_a = 20;        %Potassium (slow inactivation) parameters.
    theta_b = -70; sigma_b = -6;
    tau_a = 2; tau_b = 150;

    ERest=-70;                          %Leak reversal.
    
    tau_I = 12;                         %Time constant of inhibition.

    no_cells = size(V_in,1)/6;
    
    Vs_in = V_in(1:no_cells);
    h_in = V_in(no_cells + (1:no_cells));
    n_in = V_in(2*no_cells + (1:no_cells));
    a_in = V_in(3*no_cells + (1:no_cells));
    b_in = V_in(4*no_cells + (1:no_cells));
    s_in = V_in(5*no_cells + (1:no_cells));

    Vs_out = g_struct.gNa.*(steady(Vs_in,g_struct.theta_m,sigma_m).^3).*h_in.*(ENa-Vs_in) + g_struct.gK.*(n_in.^2).*(EK-Vs_in) ...
      + g_struct.gD.*a_in.^3.*b_in.*(EK-Vs_in) + g_struct.gL.*(ERest-Vs_in) + I_app.*(t > I_on) + noise_multiplier.*randn(no_cells,1) ...
      + CS*s_in.*(-80-Vs_in) + (CG*diag(Vs_in)-diag(Vs_in)*CG)*ones(no_cells,1);                            %Update I-cell voltage.
    h_out = (steady(Vs_in,theta_h,sigma_h)-h_in)./tau_h(Vs_in);                                             %Update h.
    n_out = (steady(Vs_in,theta_n,sigma_n)-n_in)./tau_n(Vs_in);                                             %Update n.
    a_out = (steady(Vs_in,theta_a,sigma_a)-a_in)/tau_a;                                                     %Update a.
    b_out = (steady(Vs_in,theta_b,sigma_b)-b_in)/tau_b;                                                     %Update b.
  
    s_out = ((1+tanh(Vs_in/10))/2).*(1-s_in)/0.5 - s_in/tau_I;                                              %Update s.
    
    V_out = [Vs_out; h_out; n_out; a_out; b_out; s_out];
    
end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function iout = steady(V,theta,sigma)
iout = 1./(1+exp((theta-V)/sigma));
end

function th = tau_h(V)
theta_th = -60; sigma_th = -12;
th = .5 + 14./(1+exp((theta_th-V)/sigma_th));
end

function tn = tau_n(V)
tn_1 = .087 + (11.4)./(1+exp((V+14.6)/8.6));
tn_2 = .087 + (11.4)./(1+exp((1.3-V)/18.7));
tn = tn_1.*tn_2;
end