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

function [Vs,Vd,h,n,a,b,s,t] = Golomb_2007_RK45_w_dendrite(no_cells,I_app,T0,in_g_struct,noise_multiplier,gG_multiplier,gS_multiplier)

    % Setting different conductances for soma & dendrites.
    g_mod_vec = [2*ones(no_cells,1); ones(no_cells,1)/10];

    % Default conductances.
    g_struct.gNa = 112.5*g_mod_vec;
    g_struct.gK = 225*g_mod_vec;
    g_struct.gL = 0.25*g_mod_vec;
    g_struct.g_sd = (0.25/2)*ones(no_cells,1);
    g_struct.theta_m = -24;
    g_struct.gD = 1.6;
    g_fields = fieldnames(g_struct);

    % Assigning conductances from in_g_struct.
    in_fields = fieldnames(in_g_struct);
    for field = in_fields'
        if ~any(strcmpi(char(field), g_fields))
            error(['Bad option field name: ',char(field)]);
        end
        if ~isempty(eval(['in_g_struct.',char(field)])) % empty values are discarded
            eval(['g_struct.', char(field), ' = in_g_struct.', char(field),';']);
        end
    end
    
    if isscalar(I_app)
        
        I_app = [I_app*ones(no_cells,1); zeros(no_cells,1)];
        
    end


    dt = 0.005;                       %The time step.
    T  = ceil(T0/dt);
    t = (1:T)*dt;                     %Define time axis vector (useful for plotting).

    %% Initial conditions.
    
    Vs_0 = -70+7*rand(no_cells,1);  	%Set the initial conditions for P-cells, I-cells, and synapses.
    Vd_0 = Vs_0;
    h_0 = 0.0 + 0.1*rand(2*no_cells,1);
    n_0 = 0.0 + 0.1*rand(2*no_cells,1);
    a_0 = 0.0 + 0.1*rand(2*no_cells,1);
    b_0 = 0.0 + 0.1*rand(2*no_cells,1);
    s_0 = 0.0 + 0.1*rand(2*no_cells,1);
    I_on = 100;%7*rand(no_cells,1);

    V0 = [Vs_0; Vd_0; h_0; n_0; a_0; b_0; s_0];
    
    %% Connectivity matrices.
    
    [CS,CG] = striatal_connectivity_matrices(no_cells,0,1);
    
    if nargin >= 6    % Doubling gap junction conductances, unless otherwise specified.
        
        if ~isempty(gG_multiplier)
            
            CG = gG_multiplier*CG;
            
        else
            
            CG = 2*CG;
            
        end
        
    elseif nargin < 6
        
        CG = 2*CG;
        
    end
    
    if nargin > 6 && ~isempty(gS_multiplier)
        
        CS = gS_multiplier*CS;
        
    end

    CG = [zeros(no_cells) diag(g_struct.g_sd); diag(g_struct.g_sd) CG];
    CS = [CS zeros(no_cells); zeros(no_cells,2*no_cells)];
    
    %% Perform integration using Runge-Kutta 4/5.
    [t, V] = ode45(@(t, V) Golomb_2007_RHS(t,V,CS,CG,g_struct,I_app,I_on,noise_multiplier),t,V0);%dt*[1 T],V0);
    
    Vs = V(:,1:no_cells)';
    Vd = V(:,no_cells + (1:no_cells))';
    h = V(:,2*no_cells + (1:2*no_cells))';
    n = V(:,4*no_cells + (1:2*no_cells))';
    a = V(:,6*no_cells + (1:2*no_cells))';
    b = V(:,8*no_cells + (1:2*no_cells))';
    s = V(:,10*no_cells + (1:2*no_cells))';
  
end