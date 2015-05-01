%Implement a model of a fast-spiking interneuron, as in Golomb 2007.

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

function [Vs,h,n,a,b,t] = Golomb_2007(no_cells,I0,T0,theta_m,gD,noise_multiplier)

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
    
  t = (1:T)*dt;                     %Define time axis vector (useful for plotting).
  
  Vs = zeros(no_cells,T);                %Make empty variables to hold I-cell results.
  h = zeros(no_cells,T);
  n = zeros(no_cells,T);
  a = zeros(no_cells,T);
  b = zeros(no_cells,T);
  
  Vs(:,1)=-70;%+7*rand(no_cells,1);  	%Set the initial conditions for P-cells, I-cells, and synapses.
  h(:,1)=0.0 + 0.1*rand(no_cells,1);
  n(:,1)=0.0 + 0.1*rand(no_cells,1);
  a(:,1)=0.0 + 0.1*rand(no_cells,1);
  b(:,1)=0.0 + 0.1*rand(no_cells,1);
  I_on=100;%7*rand(no_cells,1);
  
  for i=1:T-1                       %Integrate the equations.
      
      %I-cell dynamics - NOTE the synaptic current!
      Vs(:,i+1) = Vs(:,i) + dt*(gNa*(inf(Vs(:,i),theta_m,sigma_m).^3).*h(:,i).*(ENa-Vs(:,i)) + gK*(n(:,i).^2).*(EK-Vs(:,i)) ...
          + gD*a(:,i).^3.*b(:,i).*(EK-Vs(:,i)) + gL*(ERest-Vs(:,i)) + I0*(t(i)>I_on) + noise_multiplier*randn);            %Update I-cell voltage.
      h(:,i+1) = h(:,i) + dt*((inf(Vs(:,i),theta_h,sigma_h)-h(:,i))/tau_h(Vs(:,i)));                                         %Update h.
      n(:,i+1) = n(:,i) + dt*((inf(Vs(:,i),theta_n,sigma_n)-n(:,i))/tau_n(Vs(:,i)));                                         %Update n.
      a(:,i+1) = a(:,i) + dt*((inf(Vs(:,i),theta_a,sigma_a)-a(:,i))/tau_a);                                                  %Update a.
      b(:,i+1) = b(:,i) + dt*((inf(Vs(:,i),theta_b,sigma_b)-b(:,i))/tau_b);                                                  %Update b.
      
  end
  
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