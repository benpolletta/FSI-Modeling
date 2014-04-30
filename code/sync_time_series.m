function [sync_index,sync_index_times] = sync_time_series(Vs,dt)

[no_cells,length_trace]=size(Vs);

t=0:dt:(length_trace*dt-dt);

[cells,spike_indices] = find(diff(sign(diff(Vs,1,2)),1,2)<-1);

spike_times=cell(no_cells,1);
no_spikes=0;

for n=1:no_cells
    
    spike_times{n}=t(spike_indices(cells==n));
    
    no_spikes=max(no_spikes,length(spike_times{n}));
    
end

for n=1:no_cells
    
    spike_times{n}((end+1):no_spikes)=nan;
    
end

sync_index=zeros(nchoosek(no_cells,2),no_spikes);
sync_index_times=zeros(nchoosek(no_cells,2),no_spikes);

pairs=nchoosek(1:no_cells,2);

for p=1:size(pairs,1)
    
    spike_times_1=spike_times{pairs(p,1)};
    spike_times_2=spike_times{pairs(p,2)};
    
    sync_index_times(p,:)=nanmean([spike_times_1; spike_times_2]);
    
    ISI_intra=nanmean([diff(spike_times_1); diff(spike_times_2)]);
    
    ISI_inter=nansum([spike_times_1; -spike_times_2]);
    ISI_inter((length(ISI_intra)+1):end)=[];
    
    sync_index_temp=mod(ISI_inter./ISI_intra,1);
    sync_index_temp((end+1):no_spikes)=nan;
    
    sync_index(p,:)=sync_index_temp;
    
end



