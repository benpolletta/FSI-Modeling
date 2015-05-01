%Implement a model of Interneuron Network Gamma (ING).
%Use the Hodgkin-Huxley equation for the individual neuron.
%Use Borgers et al, PNAS, 2008 to model the synapses
%  (see http://www.pnas.org/content/105/46/18023.abstract).

%INPUTS:
%  no_cells = number of simulated FSIs.
%  I0  = Injected current to cell.
%  tauI = decay time of inhibitory synapse.
%  T0  = total time of simulation [s].
%  CS = inhibitory synapse connectivity matrix, including conductance (strength) for each synapse.
%  CG = gap junction connectivity matrix, including conductance (strength) for each gap junction.

%OUTPUTS:
%  V = voltage of I-cell.
%  s = inhibitory synapse.
%  t = time axis vector (useful for plotting).

function [V,s,t] = ing(no_cells,I0,tauI,T0,CS,CG)

  dt = 0.005;                       %The time step.
  T  = ceil(T0/dt);
  gNa = 120;  ENa=115;              %Sodium max conductance and reversal.    
  gK = 36;  EK=-12;                 %Potassium max conductance and reversal.
  gL=0.3;  ERest=10.6;              %Leak max conductance and reversal

  t = (1:T)*dt;                     %Define time axis vector (useful for plotting).
  
  V = zeros(no_cells,T);                %Make empty variables to hold I-cell results.
  m = zeros(no_cells,T);
  h = zeros(no_cells,T);
  n = zeros(no_cells,T);

  s = zeros(no_cells,T);                %Make empty variables to hold the synapse results.
  
  VI(:,1)=-70.0 + 10*rand(no_cells,1);  	%Set the initial conditions for P-cells, I-cells, and synapses.
  mI(:,1)=0.05 + 0.1*rand(no_cells,1);
  hI(:,1)=0.54 + 0.1*rand(no_cells,1);
  nI(:,1)=0.34 + 0.1*rand(no_cells,1);
  sI(:,1)=0.0 + 0.1*rand(no_cells,1);
  
  for i=1:T-1                       %Integrate the equations.
      
      %I-cell dynamics - NOTE the synaptic current!
      V(:,i+1) = V(:,i) + dt*(gNa*m(:,i)^3*h(:,i)*(ENa-(V(:,i)+65)) + gK*n(:,i)^4*(EK-(V(:,i)+65)) + gL*(ERest-(V(:,i)+65)) + I0 ...
          + gI*CS*s(:,i)*(-80-V(:,i))) + gG*(diag(V(:,i))*CG-CG*diag(V(:,i)))*ones(no_cells,1);                     %Update I-cell voltage.
      m(:,i+1) = m(:,i) + dt*(alphaM(V(:,i))*(1-m(:,i)) - betaM(V(:,i))*m(:,i));                                    %Update m.
      h(:,i+1) = h(:,i) + dt*(alphaH(V(:,i))*(1-h(:,i)) - betaH(V(:,i))*h(:,i));                                    %Update h.
      n(:,i+1) = n(:,i) + dt*(alphaN(V(:,i))*(1-n(:,i)) - betaN(V(:,i))*n(:,i));                                    %Update n.
      s(:,i+1) = s(:,i) + dt*(((1+tanh(V(:,i)/10))/2)*(1-s(:,i))/0.5 - s(:,i)/tauI);                                %Update s.
      
  end
  
end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function aM = alphaM(V)
aM = (2.5-0.1*(V+65)) ./ (exp(2.5-0.1*(V+65)) -1);
end

function bM = betaM(V)
bM = 4*exp(-(V+65)/18);
end

function aH = alphaH(V)
aH = 0.07*exp(-(V+65)/20);
end

function bH = betaH(V)
bH = 1./(exp(3.0-0.1*(V+65))+1);
end

function aN = alphaN(V)
aN = (0.1-0.01*(V+65)) ./ (exp(1-0.1*(V+65)) -1);
end

function bN = betaN(V)
bN = 0.125*exp(-(V+65)/80);
end
