function [data, name] = dave_RS(tspan, save_flag, varargin)

% Set tau_fast = 7, look at I_app = 2.5, ..., 3.5 to see transition from
% subthreshold oscillations to intermittent spiking to continuous spiking.

vary_label = ''; vary_cell = cell(floor(length(varargin)/3), 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/2)
        
        if isscalar(varargin{2*a})
            
            vary_label = [vary_label, sprintf('_%s_%g', varargin{2*a - 1}, varargin{2*a})];
            
        else
        
            vary_label = [vary_label, sprintf('_%s_%gto%g', varargin{2*a - 1}, varargin{2*a}(1), varargin{2*a}(end))];
        
        end
        
        vary_cell(a, :) = {'RS', varargin{2*a - 1}, varargin{2*a}};
        
    end
    
end

name = ['dave_RS', vary_label];

if size(vary_cell, 1) > 2

    no_figures = prod(cellfun(@length, vary_cell(3:end, 3)));

else
    
    no_figures = 1;
    
end

Iconst = 0; I_app = 0;


spec.populations(1).name = 'RS';
spec.populations(1).size = 1;
spec.populations(1).equations = {['V''=(I(t)+current)/Cm; V(0)=-65;',...
    'I(t)=I_app*((t/ton)*(t<=ton)+(ton<t&t<toff))*(1+rand*.25); ton=500;',... 
    'Cm_factor=.95/.25; gM_denom=.084/12.5; gM=Cm_factor*gM_denom*gNaF;',...
    'gNaP_denom=3.36; gNaP=gM/gNaP_denom;',...
    sprintf('toff=%f; I_app=%f; Iconst=%f;', tspan, I_app, Iconst),... 
    'monitor functions']};
spec.populations(1).mechanism_list = {'IBnoise','IBiNaF','IBiKDR','IBiMMichTauConst','IBleaknoisy','iNaP'};
spec.populations(1).parameters = {...
    'VIC',-65,'ICnoise',0,'Cm',.9,'El',-67,'El_std',0,'gl',.2,'Vnoise',.3,...
    'gNaF',100,'gKDR',33,'EKDR',-85,'gM',0.67,'E_M',-85,'gNaP',.2,...
    };

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