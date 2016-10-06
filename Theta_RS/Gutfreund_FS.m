function [data, name] = Gutfreund_FS(tspan, save_flag, varargin)

% Set tau_fast = 7, look at I_app = 2.5, ..., 3.5 to see transition from
% subthreshold oscillations to intermittent spiking to continuous spiking.

vary_label = ''; vary_cell = cell(length(varargin)/3, 3);

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

name = ['gutfreund_FS', vary_label];

if size(vary_cell, 1) > 2

    no_figures = prod(cellfun(@length, vary_cell(3:end, 3)));

else
    
    no_figures = 1;
    
end

I_const = 0; I_app = 0;

model_eqns = ['dv/dt=(I_const+I(t)+@current)/Cm; Cm=.25; v(0)=-65;',...
    'tau_fast=5; tau_h=tau_fast; tau_m=tau_fast;',...
    'fast_denom=1; gKDR=5/fast_denom; gNa=12.5/fast_denom;',...
    'I(t)=I_app*((t/500)*(t<=ton)+(ton<t&t<toff))*(1+rand*.25);',... % *((1-pulse/2)+pulse*(mod(t,750)<250&t>2*ton));',...
    sprintf('PPstim = 0; ton=500; toff=%f; I_app=%f; I_const=%f;', tspan, I_app, I_const),... %  (ton<t&t<toff)
    'monitor functions'];

% Setting up structure containing info for two populations.
s = [];
% RS cells.
s.populations(1).name = 'RS';
s.populations(1).size = 1;
s.populations(1).equations = model_eqns;
s.populations(1).mechanism_list = {'iNaP','iKs','iKDRG','iNaG','gleak','CaDynT','iCaT','iKCaT','iPeriodicPulses'};
%s.populations(1).parameters = {'Iapp',5,'gNa',120,'gK',36,'noise',40};
% FS cells.
s.populations(2).name = 'FS';
s.populations(2).size = 1;
s.populations(2).equations = model_eqns;
s.populations(2).mechanism_list = {'iNa','iK'};
s.populations(2).parameters = {'Iapp',0,'gNa',12.5,'gK',5};
% FS->RS GABA synapses.
s.connections(1).direction = 'FS->RS';
s.connections(1).mechanism_list = {'iGABAa'};
s.connections(1).parameters = {'tauD',10,'gSYN',.1,'netcon','ones(N_pre,N_post)'};
% RS->FS AMPA synapses.
s.connections(2).direction = 'RS->FS';
s.connections(2).mechanism_list = {'iAMPA'};
s.connections(2).parameters = {'tauD',2,'gSYN',.1,'netcon','ones(N_pre,N_post)'};

if ~isempty(varargin)
    
    % if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(s, 'tspan', [0 tspan], 'vary', vary_cell, 'parallel_flag', 1, 'downsample_factor', 25);
    
    % else
    % 
    %     data = SimulateModel(model_eqns, 'tspan', [0 tspan], 'vary', vary_cell, 'compile_flag', 1);
    % 
    % end
    
else
    
    % if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(s, 'tspan', [0 tspan], 'parallel_flag', 1, 'downsample_factor', 25);
    
    % else
    % 
    %     data=SimulateModel(model_eqns, 'tspan', [0 tspan], 'compile_flag', 1);
    % 
    % end
     
end
    
try 
    
    PlotData(data, 'variable', {'RS_v', 'FS_v'})

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