%Implement a model of Interneuron Network Gamma (ING).
%Use the Hodgkin-Huxley equation for the individual neuron.
%Use Borgers et al, PNAS, 2008 to model the synapses
%  (see http://www.pnas.org/content/105/46/18023.abstract).

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

function [Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(no_cells,I0,T0,g_sd,CS,CG)

  dt = 0.005;                       %The time step.
  T  = ceil(T0/dt);
  tauI = 12;                        %Decay time of inhibition.
  gNa = 120;  ENa=115;              %Sodium max conductance and reversal.    
  gK = 36;  EK=-12;                 %Potassium max conductance and reversal.
  gL = 0.3;  ERest=10.6;            %Leak max conductance and reversal.
  if isempty(g_sd)
    g_sd = gL/2;                    %Conductance from dendrite to soma (half of leak conductance, as in Lewis & Rinzel 2003).
  end
    
  t = (1:T)*dt;                     %Define time axis vector (useful for plotting).
  
  Vs = zeros(no_cells,T);                %Make empty variables to hold I-cell results.
  Vd = zeros(no_cells,T);
  h = zeros(no_cells,T);
  n = zeros(no_cells,T);
  a = zeros(no_cells,T);
  b = zeros(no_cells,T);

  s = zeros(no_cells,T);                %Make empty variables to hold the synapse results.
  
  Vs(:,1)=-70+70*rand(no_cells,1);  	%Set the initial conditions for P-cells, I-cells, and synapses.
  Vd(:,1)=Vs(:,1);
  h(:,1)=0.5*rand(no_cells,1);
  n(:,1)=0.35 + .4*rand(no_cells,1);
  s(:,1)=0.0 + 0.1*rand(no_cells,1);
  I_on=7*rand(no_cells,1);
  
  for i=1:T-1                       %Integrate the equations.
      
      %I-cell dynamics - NOTE the synaptic current!
      Vs(:,i+1) = Vs(:,i) + dt*(gNa*(inf(Vs(:,i),theta_m,sigma_m).^3).*h(:,i).*(ENa-(Vs(:,i)+65)) + gK*(n(:,i).^2).*(EK-(Vs(:,i)+65)) ...
          + gD*a(:,i).^3.*b(:,i).*(EK-(Vs(:,i)+65)) + gL*(ERest-(Vs(:,i)+65)) + I0*(t(i)>I_on) ...
          + CS*s(:,i).*(-80-Vs(:,i)) + g_sd*(Vd(:,i)-Vs(:,i)));                                                                                %Update I-cell voltage of soma.
      Vd(:,i+1) = Vd(:,i) + dt*(g_sd*(Vs(:,i)-Vd(:,i)) + (CG*diag(Vd(:,i))-diag(Vd(:,i))*CG)*ones(no_cells,1));
      h(:,i+1) = h(:,i) + dt*((inf(Vs(:,i),theta_h,sigma_h)-h(:,i))/tau_h(Vs(:,i)));                                                  %Update h.
      n(:,i+1) = n(:,i) + dt*((inf(Vs(:,i),theta_n,sigma_n)-n(:,i))/tau_n(Vs(:,i)));                                                  %Update n.
      a(:,i+1) = a(:,i) + dt*((inf(Vs(:,i),theta_a,sigma_a)-a(:,i))/tau_a);                                                  %Update a.
      b(:,i+1) = b(:,i) + dt*((inf(Vs(:,i),theta_b,sigma_b)-b(:,i))/tau_b);                                                  %Update b.
      
      s(:,i+1) = s(:,i) + dt*(((1+tanh(Vs(:,i)/10))/2).*(1-s(:,i))/0.5 - s(:,i)/tauI);                                                %Update s.
      
  end
  
end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function iout = inf(V,theta,sigma)
iout = 1./(1+exp((theta-V)/sigma));
end

function th = tau_h(V)
th = .5 + 14./(1+exp((theta_th-V)/sigma_th));
end

function tn = tau_n(V)
tn_1 = .087 + (11.4)./(1+exp((V+14.6)/8.6));
tn_2 = .087 + (11.4)./(1+exp((1.3-V)/18.7));
tn = tn_1.*tn_2;
end