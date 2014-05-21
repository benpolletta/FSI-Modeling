function Spike_Times = poiss_spike_times(T0,lambda,no_inputs)

% T0 should be in seconds; lambda should be in Hz.

Spike_Times = zeros(no_inputs,1);

while any(Spike_Times(:,end) < T0)
    
    ISIs = exprnd(lambda,no_inputs,round(T0*lambda));
    New_Spike_Times = cumsum([Spike_Times(:,end) ISIs],2);
    New_Spike_Times = New_Spike_Times(:,2:end);
    
    Spike_Times = [Spike_Times New_Spike_Times];
    
end

Spike_Times = Spike_Times(:,2:end);