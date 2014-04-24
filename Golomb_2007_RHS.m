function [V_out] = Golomb_2007_RHS(t,V_in,no_cells,theta_m,gD,I_app,I_on,noise_multiplier)

    gNa = 112.5;  ENa=50;             %Sodium max conductance and reversal.
    sigma_m = 11.5;
    theta_h = -58.3; sigma_h = -6.7;

    gK = 225;  EK=-90;                %Potassium max conductance and reversal.
    theta_n = -12.4; sigma_n = 6.8;   %Delayed rectifier parameters.

    theta_a = -50; sigma_a = 20;      %Potassium (slow inactivation) parameters.
    theta_b = -70; sigma_b = -6;
    tau_a = 2; tau_b = 150;

    gL = 0.25;  ERest=-70;            %Leak max conductance and reversal.

    Vs_in = V_in(1:no_cells);
    h_in = V_in(no_cells + (1:no_cells));
    n_in = V_in(2*no_cells + (1:no_cells));
    a_in = V_in(3*no_cells + (1:no_cells));
    b_in = V_in(4*no_cells + (1:no_cells));

    Vs_out = gNa*(inf(Vs_in,theta_m,sigma_m).^3).*h_in.*(ENa-Vs_in) + gK*(n_in.^2).*(EK-Vs_in) ...
      + gD.*a_in.^3.*b_in.*(EK-Vs_in) + gL*(ERest-Vs_in) + I_app*(t > I_on) + noise_multiplier*randn; %real(sqrt(2*noise_multiplier*randn));    %Update I-cell voltage.
    h_out = (inf(Vs_in,theta_h,sigma_h)-h_in)./tau_h(Vs_in);                                            %Update h.
    n_out = (inf(Vs_in,theta_n,sigma_n)-n_in)./tau_n(Vs_in);                                            %Update n.
    a_out = (inf(Vs_in,theta_a,sigma_a)-a_in)/tau_a;                                                    %Update a.
    b_out = (inf(Vs_in,theta_b,sigma_b)-b_in)/tau_b;                                                    %Update b.
  
    V_out = [Vs_out; h_out; n_out; a_out; b_out];
    
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