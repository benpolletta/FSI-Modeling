% Implement a network of fast-spiking interneurons, as in Golomb 2007, with
% connectivity as given in Russo & Taverna 2012.

%INPUTS:
%  no_cells = number of simulated FSIs.
%  I0  = Injected current to cell.
%  tauI= decay time of inhibitory synapse.
%  T0  = total time of simulation [s].
%  CS = inhibitory synapse connectivity matrix, including conductance (strength) for each synapse.
%  CG = gap junction connectivity matrix, including conductance (strength) for each gap junction.

%OUTPUTS:
%  V = voltage of I-cell.
%  s = inhibitory synapse.
%  t = time axis vector (useful for plotting).

function [Vs,Vd,h,n,a,b,t,s] = Golomb_2007_w_dendritic_gap_jxn(no_cells,I0,T0,theta_m,gD)

  dt = 0.005;                       %The time step.
  T  = ceil(T0/dt);
  
  gNa = 112.5;  ENa=50;             %Sodium max conductance and reversal.
  sigma_m = 11.5;
  theta_h = -58.3; sigma_h = -6.7;
  
  gK = 225;  EK=-90;                %Potassium max conductance and reversal.
  theta_n = -12.4; sigma_n = 6.8;   %Delayed rectifier parameters.
  
  theta_a = -50; sigma_a = 20;      %Potassium (slow inactivation) parameters.
  theta_b = -70; sigma_b = -6;
  tau_a = 2; tau_b = 150;
  
  gL = 0.25;  ERest=-70;            %Leak max conductance and reversal.
    
  g_sd = gL/2;                       %Conductance between soma and dendrite.
  
  t = (1:T)*dt;                     %Define time axis vector (useful for plotting).
  
  Vs = zeros(no_cells,T);                %Make empty variables to hold I-cell results.
  Vd = zeros(no_cells,T);
  h = zeros(no_cells,T);
  n = zeros(no_cells,T);
  a = zeros(no_cells,T);
  b = zeros(no_cells,T);
  s = zeros(no_cells,T);
  
  Vs(:,1)=-70+7*rand(no_cells,1);  	%Set the initial conditions for I-cells, and synapses.
  Vd(:,1)=Vs(:,1);
  h(:,1)=0.0 + 0.1*rand(no_cells,1);
  n(:,1)=0.0 + 0.1*rand(no_cells,1);
  a(:,1)=0.0 + 0.1*rand(no_cells,1);
  b(:,1)=0.0 + 0.1*rand(no_cells,1);
  s(:,1)=0.0 + 0.1*rand(no_cells,1);
  I_on=100;%7*rand(no_cells,1);
  
  [CS,CG] = striatal_connectivity_matrices(no_cells);
  
  for i=1:T-1      %Integrate the equations.

      Vs(:,i+1) = Vs(:,i) + dt*(gNa*(steady(Vs(:,i),theta_m,sigma_m).^3).*h(:,i).*(ENa-Vs(:,i)) + gK*(n(:,i).^2).*(EK-Vs(:,i)) ...
          + gD*a(:,i).^3.*b(:,i).*(EK-(Vs(:,i)+65)) + gL*(ERest-Vs(:,i)) + I0*(t(i)>I_on) ...
          + CS*s(:,i).*(-80-Vs(:,i)) + g_sd*(Vd(:,i)-Vs(:,i)));                                                             %Update I-cell voltage of soma.
      Vd(:,i+1) = Vd(:,i) + dt*(g_sd*(Vs(:,i)-Vd(:,i)) + (CG*diag(Vd(:,i))-diag(Vd(:,i))*CG)*ones(no_cells,1));             %Update I-cell voltage of dendrite.
      h(:,i+1) = h(:,i) + dt*((steady(Vs(:,i),theta_h,sigma_h)-h(:,i))./tau_h(Vs(:,i)));                                        %Update h.
      n(:,i+1) = n(:,i) + dt*((steady(Vs(:,i),theta_n,sigma_n)-n(:,i))./tau_n(Vs(:,i)));                                        %Update n.
      a(:,i+1) = a(:,i) + dt*((steady(Vs(:,i),theta_a,sigma_a)-a(:,i))./tau_a);                                                 %Update a.
      b(:,i+1) = b(:,i) + dt*((steady(Vs(:,i),theta_b,sigma_b)-b(:,i))./tau_b);                                                 %Update b.
      
  end
  
end

%Below, define the auxiliary functions for each gating variable.

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