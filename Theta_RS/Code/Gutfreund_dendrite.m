function data = Gutfreund_dendrite(gcom, I_app, varargin)

vary_label = ''; vary_cell = cell(length(varargin)/3, 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/3)
        
        vary_label = [vary_label, sprintf('_%s%s_%fto%f', varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a}(1), varargin{3*a}(end))];
        
        vary_cell(a, :) = {varargin{3*a - 2}, varargin{3*a - 1}, varargin{3*a}};
        
    end
    
end

eqns = ['dv/dt=I_app*(1+rand(1,N_pop)*.25)+@current/Cm; Cm=.25; v(0)=-65; I_app=0;',...
    'monitor functions'];
    % sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=500; toff=3500; I_app=%f;', I_app),...

d_num = 1; s_num = 3;

s=[];
s.populations(1).name='soma';
s.populations(1).size=1;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNaG','iKDRG','gleak','iNaP','iKs'};
s.populations(1).parameters={'gKDR',5,'gNa',12.5,'gl',0.025,'Noffset',-17.5,'Koffset',-17.5}; % ,... 'gNaP',0.025/s_num,'gKs',0.084/s_num,'halfNaP',-60,'halfKs',-50};
s.populations(2).name='dendrite';
s.populations(2).size=1;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNaG','iKDRG','gleak','iNaP','iKs'};
s.populations(2).parameters={'gNaP',0.025,'gKs',0.084,...
    'gdenomK',d_num,'gdenomN',d_num,'gdenoml',d_num}; % ,'tau_h',7,'tau_m',7,'Noffset',0,'Koffset',0};
s.connections(1).direction='soma->dendrite';
s.connections(1).mechanism_list={'iCOM'};
s.connections(1).parameters={'gCOM',gcom};
s.connections(2).direction='dendrite->soma';
s.connections(2).mechanism_list={'iCOM'};
s.connections(2).parameters={'gCOM',gcom};

if ~isempty(varargin)
    
    if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(s, 'tspan', [0 2000], 'vary', vary_cell);
    
    else
        
        data = SimulateModel(s, 'tspan', [0 2000], 'vary', vary_cell, 'compile_flag', 1);
    
    end
        
else
    
    if strcmp(version('-release'), '2012a')
    
        data = SimulateModel(s, 'tspan', [0 1500]);
    
    else
    
        data=SimulateModel(s, 'tspan', [0 1500], 'compile_flag', 1); %{'pop1','gKs',.04:.002:.06 % ;'pop1','gNaP',.010:.001:.020}); % {'pop1','gd',5:10;'pop1','I_app',10:20});

    end
        
end
    
try 
    
    PlotData(data)

    save_as_pdf(gcf, [sprintf('Figures/gutfreund_dendrite_gcom%f_Iapp%f', gcom, I_app), vary_label])

catch error
    
    display('PlotData failed:')
    
    display(error)

end