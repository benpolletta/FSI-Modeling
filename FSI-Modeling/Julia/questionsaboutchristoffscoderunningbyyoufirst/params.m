%-----------------------------------
% rng seed
rand('twister',5489); randn('state',0);  %why do we need to do this? controlling some aspect of the rng?
%-----------------------------------



%-----------------------------------
% number of E- and I-cells:
num_e=0; num_i=200;
%-----------------------------------



%-----------------------------------
% density of synaptic connections: 
p_ee=1.0; p_ei=1.0; p_ie=0.5; p_ii=0.5;
%-----------------------------------



%-----------------------------------
% density of gap-junctional connections:
p_hat_gap=0.02;
%-----------------------------------



%-----------------------------------
% rise and decay time constants associated with synapses:
tau_r_e=0.1*ones(num_e,1); tau_d_e=3*ones(num_e,1);
tau_r_i=0.3*ones(num_i,1); tau_d_i=9*ones(num_i,1); 
%-----------------------------------



%-----------------------------------
% synaptic reversal potentials: 
v_rev_e=0; v_rev_i=-80;          %i don't quite understand these
%-----------------------------------



%-----------------------------------
% time simulated (in ms):
t_final=1000;  
%-----------------------------------



%-----------------------------------
% time step used in the midpoint method:
dt=0.01; 
%-----------------------------------



%-----------------------------------
% strength of synaptic connections: 
g_hat_ie=0.0; g_hat_ei=0.0; g_hat_ii=0.5; g_hat_ee=0; %is this in nanosiemens
%-----------------------------------



%-----------------------------------
% strength of gap-junctional connections:
g_hat_gap=0.05; 					%what about this one
%-----------------------------------



%-----------------------------------
% strength and decay time constants of adaptation current:
u_g=rand(num_i,1); u_at_all=rand(num_i,1); %two random vectors, size of # of inhib cells
g_M=0.1+0.2*u_g; %random initial conditions for g_M
tau_M=100*ones(num_i,1); %this is a vector of 200 100s, not sure what for
%-----------------------------------



%-----------------------------------
% external drive to the E-cells:

% deterministic drive:
u_e=randn(num_e,1); %random vector, size of # of inhib cells
I_e=@(t) 0.2*ones(num_e,1); %I_e is a function of t, but always returns a num_e length vector of 0.2s

% maximum conductance, decay time, and frequency of Poisson train of excitatory input pulses:
g_stoch_e=0.0; f_stoch_e=40; tau_d_stoch_e=3; %that's a pretty high frequency, unless there's only one input synapse per cell, in which case it's low.
%what do we mean "conductance"? and why is it 0?
%-----------------------------------



%-----------------------------------
% external drive to the I-cells:

% deterministic drive:
a=0;  %a=0.1;
b=0;  %b=0.5;
base=2;
alpha_gamma=5;
T_gamma=20;
J_gamma=@(t) exp(alpha_gamma*(cos(pi*t/T_gamma).^2))-1; %J_gamma is a function of t. should plot
C=1/mean(J_gamma([0:1000]/1000*T_gamma));
I_gamma=@(t) C*J_gamma(t); %I_gamma(t) = C*J_gamma(t). also plot this

alpha_theta=3;
T_theta=120;
J_theta=@(t) exp(alpha_theta*(cos(pi*t/T_theta).^2))-1;
C=1/mean(J_theta([0:1000]/1000*T_theta));
I_theta=@(t) C*J_theta(t);


u_i=randn(num_i,1); %num_i length vector of random numbers
I_i=@(t) base*(1+0.1*u_i)+a*I_gamma(t)+b*I_theta(t); %base plus some gammas and thetas. what are all these gammas and thetas anyway
%wait is the external drive rhythmic? no wonder the output is rhythmic
%oh... not if a and b are zero

% maximum conductance, decay time, and frequency of Poisson train of excitatory input pulses:
if ijk==1, %so i guess this is the difference between the two runs, but i still don't get what g_stoch is
    g_stoch_i=0;
else
    g_stoch_i=0.03;
end;
f_stoch_i=20; tau_d_stoch_i=3; 
%-----------------------------------
