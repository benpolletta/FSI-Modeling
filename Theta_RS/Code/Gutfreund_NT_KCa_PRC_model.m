function [data, name] = Gutfreund_NT_KCa_PRC_model(tspan, save_flag, varargin)

% Set tau_fast = 7, look at deepRS I_app = 2.5, ..., 3.5 to see transition from
% subthreshold oscillations to intermittent spiking to continuous spiking.
% Since there are two populations, need to specify to which population (deepRS
% or deepFS) a given vary statement applies.

vary_label = ''; vary_cell = cell(length(varargin)/3, 3);

if ~isempty(varargin)

    for a = 1:(length(varargin)/3)

        if isscalar(varargin{3*a})

            vary_label = [vary_label, sprintf('_%s%s_%g', varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a})];

        else

            vary_label = [vary_label, sprintf('_%s%s_%gto%g', varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a}(1), varargin{3*a}(end))];

        end

        vary_cell(a, :) = {varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a}};

    end

end

name = ['gutfreund_NT_KCa_PRC', vary_label];

if size(vary_cell, 1) > 2

    no_figures = prod(cellfun(@length, vary_cell(3:end, 3)));

else

    no_figures = 1;

end

I_app = 0;

IC_V = -65;

spec = [];

spec.populations(1).name = 'deepFS';
spec.populations(1).size = 1;
spec.populations(1).equations = {['V''=(@current)/Cm; V(0)=' num2str(IC_V) '; monitor functions;']};
spec.populations(1).mechanism_list = {'IBitonic','IBnoise','FSiNaF','FSiKDR','IBleak'};
spec.populations(1).parameters = {...
  'VIC',-65,'ICnoise',0,'Cm',.25,'El',-67,'gl',0.1,...
  'stim',0.95,'onset',0,'offset',Inf,...
  'Vnoise',1,...
  'gNaF',100,'ENaF',50,...
  'gKDR',80,'EKDR',-95,...
  };

spec.populations(2).name = 'deepRS';
spec.populations(2).size = 1;
spec.populations(2).equations = {['V''=(Iconst+@current)/Cm; V(0)=' num2str(IC_V) '; monitor functions; monitor V.spikes(0);']};
spec.populations(2).mechanism_list = {'iNaP','iKs','iKDRG','iNaG','iLeak',...
    'CaDynT','iCaT','iKCaT','iPeriodicPulsesBen','itonicBen'}; %,'iFMPulses','iSpikeTriggeredPulse'}; % ,'iCarracedoEPSPs'}; % ,'VaryRandomSeed'}; % 'iPeriodicSpikes',
spec.populations(2).parameters = {...
  'Cm',2.7,'PPstim',0,'gSpike',0,'PPcenter',0,'PPnorm',0,...
  'gKs',0,'gNaP',1.4472/3.36,'gKCa',.1512,'bKCa',.02,'gCa',.54,'CAF',24/10.8,...
  'gl',.78,'gNa',10.8*12.5,'gKDR',10.8*5,...
  'Iconst',0,'tauh',5,'taum',5,...
  'ton',500,'toff',tspan(end),'Inoise',0.25,...                                         %  (ton<t&t<toff) %%% 'PPstim = 0; PPfreq = 1.5; PPwidth = floor((1000/PPfreq)/4); PPshift = 0; ap_pulse_num = 0; kernel_type = 7;',... % in ms
  'PPonset',500,... % 'dsfact',dsfact,...
  };

spec.connections(1).direction = 'deepRS->deepFS';
spec.connections(1).mechanism_list = {'IBaIBdbiSYNseed'};
spec.connections(1).parameters = {'gSYN',.025,'ESYN',0,'tauDx',2.5,...
    'tauRx',.25,'fanout',inf,'ICnoise',0,'gSYNhetero',0, ...
    };

spec.connections(2).direction = 'deepFS->deepRS';
spec.connections(2).mechanism_list = {'IBaIBdbiSYNseed'};
spec.connections(2).parameters = {'gSYN',.2,'ESYN',-95,'tauDx',50,...
    'tauRx',.25,'fanout',inf,'ICnoise',0,'gSYNhetero',0,...
    };

if ~isempty(varargin)

    % if strcmp(version('-release'), '2012a')

        data = dsSimulate(spec, 'tspan', [0 tspan], 'vary', vary_cell,...
            'parallel_flag', 1, 'verbose_flag', 0,... 'solver', 'euler',...
            'study_dir', ['Figures/', name]);

    % else
    %
    %     data = SimulateModel(model_eqns, 'tspan', [0 tspan], 'vary', vary_cell, 'compile_flag', 1);
    %
    % end

else

    % if strcmp(version('-release'), '2012a')

        data = dsSimulate(spec, 'tspan', [0 tspan],...
            'verbose_flag', 0,... 'solver', 'euler',...
            'study_dir', ['Figures/', name]);

    % else
    %
    %     data=SimulateModel(model_eqns, 'tspan', [0 tspan], 'compile_flag', 1);
    %
    % end

end

try

    dsPlot(data)

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

    disp('Plotting failed:')

    disp(error)

end

if save_flag

    save(['Figures/', name, '.mat'], 'data', '-v7.3')

end
