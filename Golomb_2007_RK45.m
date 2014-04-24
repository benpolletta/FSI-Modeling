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

function [Vs,h,n,a,b,t] = Golomb_2007_RK45(no_cells,I_app,T0,theta_m,gD,noise_multiplier)

    dt = 0.005;                       %The time step.
    T  = ceil(T0/dt);
    t = (1:T)*dt;

    % t = (1:T)*dt;                     %Define time axis vector (useful for plotting).

    Vs_0 = -70+7*rand(no_cells,1);  	%Set the initial conditions for P-cells, I-cells, and synapses.
    h_0 = 0.0 + 0.1*rand(no_cells,1);
    n_0 = 0.0 + 0.1*rand(no_cells,1);
    a_0 = 0.0 + 0.1*rand(no_cells,1);
    b_0 = 0.0 + 0.1*rand(no_cells,1);
    I_on = 100;%7*rand(no_cells,1);

    V0 = [Vs_0; h_0; n_0; a_0; b_0];

    % Perform integration using Runge-Kutta 4/5.
    [t, V] = ode45(@(t, V) Golomb_2007_RHS(t,V,no_cells,theta_m,gD,I_app,I_on,noise_multiplier),t,V0);%dt*[1 T],V0);
    
    Vs = V(:,1:no_cells)';
    h = V(:,no_cells + (1:no_cells))';
    n = V(:,2*no_cells + (1:no_cells))';
    a = V(:,3*no_cells + (1:no_cells))';
    b = V(:,4*no_cells + (1:no_cells))';
  
end