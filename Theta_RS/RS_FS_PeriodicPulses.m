function [data, name] = RS_FS_PeriodicPulses(I_const, tspan, save_flag, Nrs, Nfs, varargin)

if isempty(Nfs), Nfs = 25; end

if isempty(Nrs), Nrs = 25; end

vary_label = ''; vary_cell = cell(floor(length(varargin)/3), 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/3)
        
        if isscalar(varargin{3*a})
            
            vary_label = [vary_label, sprintf('_%s_%s_%g', varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a})];
            
        else
        
            vary_label = [vary_label, sprintf('_%s_%s_%gto%g', varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a}(1), varargin{3*a}(end))];
        
        end
        
        vary_cell(a, :) = {varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a}};
        
    end
    
end

name = ['RS_FS_periodicpulses', vary_label];

if size(vary_cell, 1) > 2

    no_figures = prod(cellfun(@length, vary_cell(3:end, 3)));

else
    
    no_figures = 1;
    
end

% I_app = 0;

C_mult = 2.7/.25;

ton = 500;

% % % % %  RS Cells  % % % % %
spec.populations(1).name = 'RS';
spec.populations(1).size = Nrs;
spec.populations(1).equations = ['dv/dt=(I_const+@current)/Cm; Cm=2.7; v(0)=-65;',... % +I(t)
    'gNaP_denom=3.36; gNaP=gKs/gNaP_denom;',...
    sprintf('I_const=%g;', I_const),...    %  halfKs=-60; halfNaP=-60; 
    'tau_fast=5; tau_h=tau_fast; tau_m=tau_fast;',... %%% 'offset=0; Koffset=offset; Noffset=offset;',... %'slow_offset=0; halfKs=-35-slow_offset; halfNaP=-40-slow_offset;',... 'fast_denom=1; gKDR=23.4/fast_denom; gNa=45/fast_denom;',... % C_mult*ones(1, 2)),...
    'monitor functions'];
spec.populations(1).mechanism_list = {'iNaP','iKs','iNaG','iKDRG','gleak','CaDynT','iCaT','iKCaT','iStep','iPeriodicPulses'};
spec.populations(1).parameters = {...
    'gNa', C_mult*12.5, 'gKDR', C_mult*5, 'gl', C_mult*.025, 'CAF', 24/C_mult,...
    'tau_h', 5, 'tau_m', 5,...
    'gKs', 1.4472, 'gCa', .864, 'gKCa', .216, 'bKCa', .002,... % 'tau_div', 2,... % C_mult*ones(1, 3)),... %%% 
    'ton', ton, 'toff', tspan,... % 'I_app', I_app,... %  (ton<t&t<toff) %%% 'PPstim = 0; PPfreq = 1.5; PPwidth = floor((1000/PPfreq)/4); PPshift = 0; ap_pulse_num = 0; kernel_type = 7;',... % in ms % 'noise=.25; I(t)=C_mult*I_app*((t/ton)*(t<=ton)+(ton<t&t<toff))*(1+rand*noise);',... % *((1-pulse/2)+pulse*(mod(t,750)<250&t>2*ton));',...
    };

% % % % %  FS Cells  % % % % %
spec.populations(2).name = 'FS';
spec.populations(2).size = Nfs;
spec.populations(2).equations = 'dv/dt=(+@current)/Cm; Cm = .9; v(0)=-65; monitor functions';
spec.populations(2).mechanism_list = {'FSiNaF','FSiKDR','gleak'};
spec.populations(2).parameters = {'gNaF', 100, 'gKDR', 80, 'gl', .2};

% % % % % % % % % % % %  Connections  % % % % % % % % % % % % %  

i=0;

% % % % % Connection Parameters % % % % %


% % % % % % % % % % % % %  Synaptic connections % % % % % % % % % % % % %  

% Gap junction connections
ggjaRS=.0/Nrs;  % RS -> RS % Disabled RS-RS gap junctions because otherwise the Eleaknoise doesn't have any effect
ggjFS=.2/Nfs;  % IBa -> IBa

% Synapse heterogenity
gsyn_hetero = 0;
    
% RS-FS circuit (deep connections)
gAMPA_rsrs=0.1/Nrs;
gNMDA_RSRS=5/Nrs;
gAMPA_rsfs=0.4/Nrs;
gNMDA_rsfs=0/Nrs;
gGABAaff=.5/Nfs;
gGABAa_fsrs=10^(-1)/Nfs;

% Synaptic time constants & reversals
tauAMPAr=.25;  % ms, AMPA rise time; Jung et al
tauAMPAd=1;   % ms, AMPA decay time; Jung et al
tauNMDAr=5; % ms, NMDA rise time; Jung et al
tauNMDAd=100; % ms, NMDA decay time; Jung et al
tauGABAar=.5;  % ms, GABAa rise time; Jung et al
tauGABAad=8;   % ms, GABAa decay time; Jung et al
% tauGABAbr=38;  % ms, GABAa rise time; From NEURON Delta simulation
% tauGABAbd=150;   % ms, GABAa decay time; From NEURON Delta simulation
EAMPA=0;
EGABA=-95;

% % RS->RS recurrent synaptic and gap connections
if Nrs > 1
    i=i+1;
    spec.connections(i).direction = 'RS->RS';
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed','iNMDA','IBaIBaiGAP'};
    spec.connections(i).parameters = {'g_SYN',gAMPA_rsrs,'E_SYN',EAMPA,'tauDx',tauAMPAd, ...
        'tauRx',tauAMPAr,'fanout',inf,'IC_noise',0,'g_SYN_hetero',gsyn_hetero, ...
        'gNMDA',gNMDA_RSRS,'ENMDA',EAMPA,'tauNMDAr',tauNMDAr,'tauNMDAd',tauNMDAd ...
        'g_GAP',ggjaRS, ...
        };
end

% % RS->FS synaptic connection
if Nrs > 0 && Nfs > 0
    i=i+1;
    spec.connections(i).direction = 'RS->FS';
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed','iNMDA'};
    spec.connections(i).parameters = {'g_SYN',gAMPA_rsfs,'E_SYN',EAMPA,'tauDx',tauAMPAd, ...
        'tauRx',tauAMPAr,'fanout',inf,'IC_noise',0,'g_SYN_hetero',gsyn_hetero, ...
        'gNMDA',gNMDA_rsfs,'ENMDA',EAMPA,'tauNMDAr',tauNMDAr,'tauNMDAd',tauNMDAd ...
        };
end

% % % % %  FS Cells  % % % % %
% % FS->FS Synaptic connections
if Nfs > 0
    i=i+1;
    spec.connections(i).direction = 'FS->FS';                   % GABA_A
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed','IBaIBaiGAP'};
    spec.connections(i).parameters = {'g_SYN',gGABAaff,'E_SYN',EGABA,'tauDx',tauGABAad,...
        'tauRx',tauGABAar,'fanout',inf,'IC_noise',0,'g_SYN_hetero',gsyn_hetero,...
        'g_GAP',ggjFS,...
        };
end

% % FS->RS Synaptic connections
if Nfs > 0 && Nrs > 0
    i=i+1;
    spec.connections(i).direction = 'FS->RS';                   % GABA_A
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed'};
    spec.connections(i).parameters = {'g_SYN',gGABAa_fsrs,'E_SYN',EGABA,...
        'tauDx',tauGABAad,'tauRx',tauGABAar,'fanout',inf,'IC_noise',0,'g_SYN_hetero',gsyn_hetero,...
        };
end

if ~isempty(varargin)
    
    % if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(spec, 'tspan', [0 tspan], 'vary', vary_cell, 'parallel_flag', 1, 'downsample_factor', 25);
    
    % else
    % 
    %     data = SimulateModel(model_eqns, 'tspan', [0 tspan], 'vary', vary_cell, 'compile_flag', 1);
    % 
    % end
    
else
    
    % if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(spec, 'tspan', [0 tspan], 'parallel_flag', 1, 'downsample_factor', 25);
    
    % else
    % 
    %     data=SimulateModel(model_eqns, 'tspan', [0 tspan], 'compile_flag', 1);
    % 
    % end
     
end
    
try 
    
    PlotData(data)

    if no_figures > 1
        
        for f = 1:no_figures
            
            save_as_pdf(f, ['Figures/', name, sprintf('_%g', f)], '-v7.3')
            
        end
        
    else
    
        save_as_pdf(gcf, ['Figures/', name], '-v7.3')

    end
    
catch error
    
    display('PlotData failed:')
    
    display(error)

end

if save_flag
    
    save([name, '.mat'], 'data', '-v7.3')
    
end