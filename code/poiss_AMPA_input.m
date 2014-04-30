function AMPA_conductances=poiss_AMPA_input(dt,T0,lambda,no_inputs)

LT=ceil(T0/dt);
t_vec=(1:LT)*dt;
EPSP=.12*exp(-2*t_vec);

AMPA_conductances=zeros(no_inputs,LT);

% for i=1:no_inputs
% 
%     spike_times=[];
%     
%     Time=0;
%     
%     while(Time<T0)
%         
%         new_ISIs=exprnd(lambda,1,lambda*T0);
%         new_spike_times=cumsum([Time new_ISIs],2);
%         spike_times=[spike_times new_spike_times(2:end)];
%         
%         Time=spike_times(end);
% 
%     end
%     
%     no_spikes=length(spike_times);
%     
%     for s=1:no_spikes
%         
%         AMPA_conductances(i,t_vec>spike_times(s))=EPSP(1:sum(t_vec>spike_times(s)));
%             
%     end
%     
% end

% while sum(T_init<T0)>0
%     
%     Spike_Times=exprnd(lambda,no_inputs,lambda*T0);
%     Spike_Times=cumsum([T_init Spike_Times],2);
%     Spike_Times=Spike_Times(:,2:end);
%     
%     T_init=Spike_Times(end);
%     
% end

for i=1:no_inputs

    Time=0;
    
    while(Time<T0)
        
        ISI=-log(rand)/lambda;
        
        if ~isinf(ISI)
    
            Time=Time+ISI;
        
            AMPA_conductances(i,t_vec>Time)=EPSP(1:sum(t_vec>Time));
            
        end
            
%         EPSP=zeros(size(t_vec));
%         EPSP(t_vec>Time)=.12*exp(-2*(t_vec(t_vec>Time)-Time));
%         
%         AMPA_conductances(i,:)=AMPA_conductances(i,:)+EPSP;

    end
    
end