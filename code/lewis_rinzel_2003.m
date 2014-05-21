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

function [Vs,Vd,m,h,n31_32,n13,s,t] = lewis_rinzel_2003(no_cells,I0,T0,CS,CG)

  dt = 0.005;                       %The time step.
  T  = ceil(T0/dt);
  tauI = 10;                        %Decay time of inhibition.
  gNa = 900;  ENa=115;              %Sodium max conductance and reversal.    
  gK31_32 = 1800;  EK=-85;          %Potassium max conductance and reversal.
  gK13 = 1.8;
  gL = 10;  ERest=-70;            %Leak max conductance and reversal.
  g_sd = 5;                      %Conductance from dendrite to soma (half of leak conductance, as in Lewis & Rinzel 2003).

  t = (1:T)*dt;                     %Define time axis vector (useful for plotting).
  
  Vs = zeros(no_cells,T);                %Make empty variables to hold I-cell results.
  Vd = zeros(no_cells,T);
  m = zeros(no_cells,T);
  h = zeros(no_cells,T);
  n31_32 = zeros(no_cells,T);
  n13 = zeros(no_cells,T);

  s = zeros(no_cells,T);                %Make empty variables to hold the synapse results.
  
  Vs(:,1)=-70+10*rand(no_cells,1);  	%Set the initial conditions for P-cells, I-cells, and synapses.
  Vd(:,1)=Vs(:,1);
%   m(:,1)=0.05 + 0.1*rand(no_cells,1);
%   h(:,1)=0.54 + 0.1*rand(no_cells,1);
%   n31_32(:,1)=0.34 + 0.1*rand(no_cells,1);
%   n13(:,1)=0.34 + 0.1*rand(no_cells,1);
%   s(:,1)=0.0 + 0.1*rand(no_cells,1);
  
  for i=1:T-1                       %Integrate the equations.
      
      %I-cell dynamics - NOTE the synaptic current!
      Vs(:,i+1) = Vs(:,i) + dt*(gNa*(m(:,i).^3).*h(:,i).*(ENa-(Vs(:,i)+65)) + gK31_32*(n31_32(:,i).^2).*(EK-(Vs(:,i)+65)) ...
          + gK13*(n13(:,i).^4).*(EK-(Vs(:,i)+65)) + gL*(ERest-(Vs(:,i)+65)) + I0 ...
          + CS*s(:,i).*(-80-Vs(:,i)) + g_sd*(Vd(:,i)-Vs(:,i)));                                                           %Update I-cell voltage of soma.
      Vd(:,i+1) = Vd(:,i) + dt*(g_sd*(Vs(:,i)-Vd(:,i)) + (CG*diag(Vd(:,i))-diag(Vd(:,i))*CG)*ones(no_cells,1));                     
    %   m(:,i+1) = m(:,i) + dt*(alphaM(Vs(:,i)).*(1-m(:,i)) - betaM(Vs(:,i)).*m(:,i));                                    %Update m.
    %   h(:,i+1) = h(:,i) + dt*(alphaH(Vs(:,i)).*(1-h(:,i)) - betaH(Vs(:,i)).*h(:,i));                                    %Update h.
    %   n31_32(:,i+1) = n31_32(:,i) + dt*(alphaN31_32(Vs(:,i)).*(1-n31_32(:,i)) - betaN31_32(Vs(:,i)).*n31_32(:,i));      %Update n.
    %   n13(:,i+1) = n13(:,i) + dt*(alphaN13(Vs(:,i)).*(1-n13(:,i)) - betaN13(Vs(:,i)).*n13(:,i));
      m(:,i+1) = m(:,i) + dt*(steady(alphaM(Vs(:,i)), betaM(Vs(:,i))) - m(:,i))./tau(alphaM(Vs(:,i)), betaM(Vs(:,i))); %Update m.
      h(:,i+1) = h(:,i) + dt*(steady(alphaH(Vs(:,i)), betaH(Vs(:,i))) - h(:,i))./tau(alphaH(Vs(:,i)), betaH(Vs(:,i))); %Update h.
      n31_32(:,i+1) = n31_32(:,i) + dt*(alphaN31_32(Vs(:,i)).*(1-n31_32(:,i)) - betaN31_32(Vs(:,i)).*n31_32(:,i));               %Update n.
      n13(:,i+1) = n13(:,i) + dt*(alphaN13(Vs(:,i)).*(1-n13(:,i)) - betaN13(Vs(:,i)).*n13(:,i));       
      s(:,i+1) = s(:,i) + dt*(((1+tanh(Vs(:,i)/10))/2).*(1-s(:,i))/0.5 - s(:,i)/tauI);                                 %Update s.
      
  end
  
end

%Below, define the auxiliary functions alpha & beta for each gating variable.

function s = steady(alpha,beta)
s = alpha./(alpha+beta);
end

function t = tau(alpha,beta)
t = 1./(alpha + beta);
end

% function aM = alphaM(V)
% aM = (3020 - 40*V)./(exp((75 + V)/(13.5)) - 1);
% end
% 
% function bM = betaM(V)
% bM = 1.2262./exp(V/42.248);
% end

function aM = alphaM(V)
aM = (2.5-0.1*(V+65)) ./ (exp(2.5-0.1*(V+65)) -1);
end

function bM = betaM(V)
bM = 4*exp(-(V+65)/18);
end

% function aH = alphaH(V)
% aH = 0.0035./exp(V/24.186);
% end
% 
% function bH = betaH(V)
% bH = -(0.8712 + 0.017*V)./(exp((51.25 + V)/(-5.2)) - 1);
% end

function aH = alphaH(V)
aH = 0.07*exp(-(V+65)/20);
end

function bH = betaH(V)
bH = 1./(exp(3.0-0.1*(V+65))+1);
end

function aN = alphaN31_32(V)
aN = (95 - V)./(exp((V - 95)/(-11.8)) - 1);
end

function bN = betaN31_32(V)
bN = 0.025./exp(V/22.222);
end

function aN = alphaN13(V)
aN = -(0.616 + 0.014*V)./(exp((44 + V)/(-2.3)) - 1);
end

function bN = betaN13(V)
bN = 0.0043./exp((44 + V)/34);
end
