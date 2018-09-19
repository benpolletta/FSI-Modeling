function [data, name] = dave_RS_FS(tspan, Nrs, Nfs, noise, synapses, save_flag, varargin)
% Model: Kramer 2008, PLoS Comp Bio

tic

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

if ~noise, noise_flag = '_no_noise'; else noise_flag = ''; end

if ~synapses, syn_flag = '_no_syn'; else syn_flag = ''; end

name = ['dave_RS_FS', noise_flag, syn_flag, vary_label];

% Simulation mode
% sim_mode = 1;   % 1 - normal sim
%                 % 2 - sim study IB disconnected; iM and iCaH
%                 % 3 - sim study IB disconnected; current injection
%                 % 4 - sim study IB connected; vary AMPA, NMDA injection
%                 % 5 - sim study IB connected; gamma input
%                 % 6 - sim study IB connected; single pulse
%                 % 7 - sim study NG; gamma input
%                 % 8 - sim study Median Nerve phase
%                 % 9 - sim study FS-RS circuit vary RS stim
%                 % 10 - Vary iPeriodicPulses in all cells
                
                
% Cells to include in model
include_RS = 1;
include_FS = 1;

% Simulation switches
if isempty(noise), noise = 1; end
if isempty(synapses), synapses = 1; end

% number of cells per population
if isempty(Nrs), Nrs=24; end % Number of RS cells
if isempty(Nfs), Nfs=24; end % Number of FS cells

% % % % % % % % % % % % %  Injected currents % % % % % % % % % % % % %  
% tonic input currents
Jfs = 1.25;     % FS current injection; step1
JRS1 = 3;
JRS2 = 1.7;

RS_offset1=0;
RS_onset2=0;

% Poisson IPSPs to IBdb (basal dendrite)
ERAN=0;
tauRAN=2;
lambda = 1000;
RSgRAN=0.005;


% % Periodic pulse stimulation
pulse_mode = 0;
switch pulse_mode
    case 0                  % No stimulation
        PPfreq = 4; % in Hz
        PPwidth = 2; % in ms
        PPshift = 0; % in ms
        PPonset = 10;    % ms, onset time
        PPoffset = tspan(end)-0;   % ms, offset time
        %PPoffset=270;   % ms, offset time
        ap_pulse_num = 0;        % The pulse number that should be delayed. 0 for no aperiodicity.
        ap_pulse_delay = 0;  % ms, the amount the spike should be delayed. 0 for no aperiodicity.
        width2_rise = 0.75;  % Not used for Gaussian pulse
        kernel_type = 2;
        RSPPstim = 0;
        FSPPstim = 0;
    case 1                  % Gamma stimulation
        PPfreq = 40; % in Hz
        PPwidth = 2; % in ms
        PPshift = 0; % in ms
        PPonset = 250;    % ms, onset time
        PPoffset = tspan(end)-250;   % ms, offset time
        %PPoffset=270;   % ms, offset time
        ap_pulse_num = 0;        % The pulse number that should be delayed. 0 for no aperiodicity.
        ap_pulse_delay = 11;  % ms, the amount the spike should be delayed. 0 for no aperiodicity.
        width2_rise = .5;  % Not used for Gaussian pulse
        kernel_type = 1;
        RSPPstim = -4;
        FSPPstim = 0;
%         FSPPstim = -5;
%         supRSPPstim = -7;

end


% % % % % % % % % % % % %  Synaptic connections % % % % % % % % % % % % %  

% Gap junction connection
% % Deep cells
ggjaRS=.0/Nrs;  % RS -> RS % Disabled RS-RS gap junctions because otherwise the Eleaknoise doesn't have any effect
ggjFS=.2/Nfs;  % IBa -> IBa

% Synapse heterogenity
gSYNhetero = 0;

if synapses
    
    % RS-FS circuit (deep connections)
    % #mysynapses
    gAMPA_rsrs=0.1/Nrs;
    gNMDA_RSRS=5/Nrs;
    gAMPA_rsfs=0.4/Nrs;
    gNMDA_rsfs=0/Nrs;
    gGABAaff=.5/Nfs;
    gGABAa_fsrs=1.0/Nfs;
    
end

% % % % % % % % % % % % % % % % % % % % % % 

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
% TmaxGABAB=0.5;



% Current injection noise levels
FS_Vnoise = 3;
RSda_Vnoise = .3;


% constant biophysical parameters
Cm=.9;        % membrane capacitance
gl=.1;
ENa=50;      % sodium reversal potential
E_EKDR=-95;  % potassium reversal potential for excitatory cells
ECa=125;     % calcium reversal potential
ICnoise=.25;% fractional noise in initial conditions
El_std = 10;

if ~noise
    
    ICnoise = 0;
    FS_Vnoise = 0;
    RSda_Vnoise = 0;
    ERAN = 0;
    RSgRAN = 0;
    El_std = 0;
    
end

IC_V = -65;


% % % % % % % % % % % % % % % % % Override some defaults % % % % % % % % % % 
% switch sim_mode
%     case 1                                                                  % Everything default, single simulation
%         
%         vary = [];
%     
%     case 9  % Vary RS cells in RS-FS network
% 
%         vary = { %'RS','stim2',linspace(3.5,1.5,4); ...
%                  %'RS','PPstim',linspace(-5,-2,4); ...
%                  'RS->FS','gSYN',[.2:.1:.5]/Nrs;...
%                  'FS->RS','gSYN',[.6:.1:1.2]/Nfs;...
% 
%                  }; 
% 
%         
% end

% % % % % % % % % % % % %  Populations  % % % % % % % % % % % % %  

spec=[];
i=0;


if include_RS
    i=i+1;
    spec.populations(i).name = 'RS';
    spec.populations(i).size = Nrs;
    spec.populations(i).equations = {['V''=(current)/Cm; V(0)=' num2str(IC_V), '; monitor functions']};
    spec.populations(i).mechanism_list = {'iPeriodicPulses','IBdbiPoissonExpJason','itonicPaired','IBnoise','IBiNaF','IBiKDR','IBiMMich','IBiCaH','IBleaknoisy'};
    spec.populations(i).parameters = {...
      'VIC',-65,'ICnoise',ICnoise,'Cm',Cm,'El',-67,'El_std',El_std,'gl',gl,...
      'PPstim', RSPPstim, 'PPfreq', PPfreq,'PPwidth',PPwidth,'PPshift',PPshift,...
      'PPonset', PPonset, 'PPoffset', PPoffset,'ap_pulse_num', ap_pulse_num,...
      'ap_pulse_delay', ap_pulse_delay,'kernel_type', kernel_type, 'width2_rise', width2_rise,...
      'gRAN',RSgRAN,'ERAN',ERAN,'tauRAN',tauRAN,'lambda',lambda,...
      'stim',JRS1,'onset',0,'offset',RS_offset1,'stim2',JRS2,'onset2',RS_onset2,'offset2',Inf,...
      'Vnoise',RSda_Vnoise,...
      'gNaF',100,'ENaF',ENa,...
      'gKDR',80,'EKDR',E_EKDR,...
      'gM',0.5,'E_M',E_EKDR,...
      'gCaH',0,'E_CaH',ECa,...
      };
end


if include_FS
    i=i+1;
    spec.populations(i).name = 'FS';
    spec.populations(i).size = Nfs;
    spec.populations(i).equations = {['V''=(current)/Cm; V(0)=' num2str(IC_V) ,'; monitor functions']};
    spec.populations(i).mechanism_list = {'iPeriodicPulses','IBitonic','IBnoise','FSiNaF','FSiKDR','IBleaknoisy'};
    spec.populations(i).parameters = {...
      'VIC',-65,'ICnoise',ICnoise,'Cm',Cm,'El',-67,'El_std',20,'gl',0.1,...
      'PPstim',FSPPstim,'PPfreq',PPfreq,'PPwidth',PPwidth,'PPshift',PPshift,...
      'PPonset',PPonset,'PPoffset',PPoffset,'ap_pulse_num',ap_pulse_num,...
      'ap_pulse_delay',ap_pulse_delay,'kernel_type', kernel_type, 'width2_rise', width2_rise,...
      'stim',Jfs,'onset',0,'offset',Inf,...
      'Vnoise',FS_Vnoise,...
      'gNaF',100,'ENaF',ENa,...
      'gKDR',80,'EKDR',E_EKDR,...
      };
end

% % % % % % % % % % % %  Connections  % % % % % % % % % % % % %  

i=0;

% % % % %  RS Cells  % % % % %
% % RS->RS recurrent synaptic and gap connections
if include_RS
    i=i+1;
    spec.connections(i).direction = 'RS->RS';
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed','iNMDA','IBaIBaiGAP'};
    spec.connections(i).parameters = {'gSYN',gAMPA_rsrs,'ESYN',EAMPA,'tauDx',tauAMPAd,'tauRx',tauAMPAr,'fanout',inf,'ICnoise',0,'gSYNhetero',gSYNhetero, ...
        'gNMDA',gNMDA_RSRS,'ENMDA',EAMPA,'tauNMDAr',tauNMDAr,'tauNMDAd',tauNMDAd ...
        'g_GAP',ggjaRS,...
        };
end

% % RS->FS synaptic connection
if include_RS && include_FS
    i=i+1;
    spec.connections(i).direction = 'RS->FS';
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed','iNMDA'};
    spec.connections(i).parameters = {'gSYN',gAMPA_rsfs,'ESYN',EAMPA,'tauDx',tauAMPAd,'tauRx',tauAMPAr,'fanout',inf,'ICnoise',0,'gSYNhetero',gSYNhetero, ...
        'gNMDA',gNMDA_rsfs,'ENMDA',EAMPA,'tauNMDAr',tauNMDAr,'tauNMDAd',tauNMDAd ...
        };
end

% % % % %  FS Cells  % % % % %
% % FS->FS Synaptic connections
if include_FS
    i=i+1;
    spec.connections(i).direction = 'FS->FS';                   % GABA_A
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed','IBaIBaiGAP'};
    spec.connections(i).parameters = {'gSYN',gGABAaff,'ESYN',EGABA,'tauDx',tauGABAad,'tauRx',tauGABAar,'fanout',inf,'ICnoise',0,'gSYNhetero',gSYNhetero,...
        'g_GAP',ggjFS,...
        };
end

% % FS->RS Synaptic connections
if include_FS && include_RS
    i=i+1;
    spec.connections(i).direction = 'FS->RS';                   % GABA_A
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed'};
    spec.connections(i).parameters = {'gSYN',gGABAa_fsrs,'ESYN',EGABA,...
        'tauDx',tauGABAad,'tauRx',tauGABAar,'fanout',inf,'ICnoise',0,'gSYNhetero',gSYNhetero,...
        };
end


% % % % % % % % % % % %  Run simulation  % % % % % % % % % % % % % 

% % simulation controls
% tspan=[0 tspan]; dt=.01; solver='euler'; % euler, rk2, rk4
% dsfact=25; % downsample factor, applied after simulation

if ~isempty(varargin)
    
    data = SimulateModel(spec, 'tspan', [0 tspan], 'vary', vary_cell, 'parallel_flag', 1, 'downsample_factor', 25, 'random_seed', 1);
    
else
    
    data = SimulateModel(spec, 'tspan', [0 tspan], 'parallel_flag', 1, 'downsample_factor', 25, 'random_seed', 1);
    
end

% 'tspan',tspan,'dt',dt,'downsample_factor',dsfact,'solver',solver,... % 'coder',0,'random_seed',1,'compile_flag',1,'vary',vary,'parallel_flag',double(sim_mode ~= 1),'verbose_flag',1);
% SimulateModel(spec,'tspan',tspan,'dt',dt,'dsfact',dsfact,'solver',solver,'coder',0,'random_seed',1,'compile_flag',1,'vary',vary,'parallel_flag',0,...
%     'cluster_flag',1,'save_data_flag',1,'study_dir','kramerout_cluster_2','verbose_flag',1);

% % Downsample data
% data = DownsampleData(data,max(round(0.1/dt),1));   % Downsample so that sampling rate is 10000 Hz (dt = 0.1 ms)


% Calculate averages across cells (e.g. mean field)
% data2 = CalcAverages(data);

% % % % % % % % % % % %  Plotting  % % % % % % % % % % % % % 
% switch sim_mode
%     case 1
%         PlotData(data,'plot_type','waveform');
% %          PlotData(data,'plot_type','rastergram');
%         
%         % PlotFR2(data);
%     
%     case 9
%         %%
%         %PlotData(data,'plot_type','waveform');
%         %PlotData(data,'plot_type','power');
%         
%         PlotData(data2,'plot_type','waveform','variable','FS_FS_IBaIBdbiSYNseed_s');
%         PlotData(data2,'plot_type','waveform','variable','RS_V');
%         PlotData(data2,'plot_type','waveform','variable','FS_V');
% 
%         PlotData(data,'plot_type','rastergram','variable','RS_V');
%         PlotData(data,'plot_type','rastergram','variable','FS_V');
%         % PlotFR2(data,'variable','RS_V'); 
%         % PlotFR2(data,'variable','FS_V'); 
%         % PlotFR2(data,'variable','RS_V','plot_type','meanFR');
%         % PlotFR2(data,'variable','FS_V','plot_type','meanFR');
%         
%     otherwise
%         PlotData(data,'plot_type','waveform');
% end

    
try 
    
    PlotData(data)

    % if no_figures > 1
    % 
    %     for f = 1:no_figures
    % 
    %         save_as_pdf(f, ['Figures/', name, sprintf('_%g', f)], '-v7.3')
    % 
    %     end
    % 
    % else
    
        save_as_pdf(gcf, ['Figures/', name], '-v7.3')

    % end
    
catch error
    
    display('PlotData failed:')
    
    display(error)

end

if save_flag
    
    save([name, '.mat'], 'data', '-v7.3')
    
end

toc

%%




