function [V_out] = Golomb_2007_RHS_w_dendrite(t,V_in,no_cells,CS,CG,theta_m,gD_in,I_app,I_on,noise_multiplier)

    gD = 2*gD_in;                       %Doubling conductances to account for dendrite.

    gNa = 2*112.5;  ENa=50;             %Sodium max conductance and reversal.
    sigma_m = 11.5;
    theta_h = -58.3; sigma_h = -6.7;

    gK = 2*225;  EK=-90;                %Potassium max conductance and reversal.
    theta_n = -12.4; sigma_n = 6.8;     %Delayed rectifier parameters.

    theta_a = -50; sigma_a = 20;        %Potassium (slow inactivation) parameters.
    theta_b = -70; sigma_b = -6;
    tau_a = 2; tau_b = 150;

    gL = 2*0.25;  ERest=-70;            %Leak max conductance and reversal.
    
    g_sd = gL/4;                        %Conductance of soma/dendrite connection.
    
    tau_I = 12;                         %Decay time of inhibition.

    Vs_in = V_in(1:no_cells);
    Vd_in = V_in(no_cells + (1:no_cells));
    h_in = V_in(2*no_cells + (1:2*no_cells));
    n_in = V_in(4*no_cells + (1:2*no_cells));
    a_in = V_in(6*no_cells + (1:2*no_cells));
    b_in = V_in(8*no_cells + (1:2*no_cells));
    s_in = V_in(10*no_cells + (1:no_cells));

    Vs_out = gNa*(inf(Vs_in,theta_m,sigma_m).^3).*h_in(1:no_cells).*(ENa-Vs_in) + gK*(n_in(1:no_cells).^2).*(EK-Vs_in) ...
      + gD.*a_in(1:no_cells).^3.*b_in(1:no_cells).*(EK-Vs_in) + gL*(ERest-Vs_in) + I_app*(t > I_on) + noise_multiplier*randn(no_cells,1) ...
      + CS*s_in.*(-80-Vs_in) + g_sd*(Vd_in-Vs_in);                                                      %Update voltage of soma.
    Vd_out = (gNa/10)*(inf(Vd_in,theta_m,sigma_m).^3).*h_in(no_cells+(1:no_cells)).*(ENa-Vd_in) + (gK/10)*(n_in(no_cells+(1:no_cells)).^2).*(EK-Vd_in) ...
      + (gD/10).*a_in(no_cells+(1:no_cells)).^3.*b_in(no_cells+(1:no_cells)).*(EK-Vd_in) + (gL/10)*(ERest-Vd_in) ...
      + g_sd*(Vs_in-Vd_in) + (CG*diag(Vd_in)-diag(Vd_in)*CG)*ones(no_cells,1);                          %Update voltage of dendrite.
    h_out = (inf(V_in(1:2*no_cells),theta_h,sigma_h)-h_in)./tau_h(V_in(1:2*no_cells));                                            %Update h.
    n_out = (inf(V_in(1:2*no_cells),theta_n,sigma_n)-n_in)./tau_n(V_in(1:2*no_cells));                                            %Update n.
    a_out = (inf(V_in(1:2*no_cells),theta_a,sigma_a)-a_in)/tau_a;                                                    %Update a.
    b_out = (inf(V_in(1:2*no_cells),theta_b,sigma_b)-b_in)/tau_b;                                                    %Update b.
  
    s_out = ((1+tanh(Vs_in/10))/2).*(1-s_in)/0.5 - s_in/tau_I;                                           %Update s.
    
    V_out = [Vs_out; Vd_out; h_out; n_out; a_out; b_out; s_out];
    
end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function iout = inf(V,theta,sigma)
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