function ISIm = mean_ISI(Vs,dt)

[no_cells,length_trace]=size(Vs);

t=0:dt:(length_trace*dt-dt);

[cells,spike_indices] = find(diff(sign(diff(Vs,1,2)),1,2)<-1);

ISIm=zeros(1,nchoosek(no_cells,2));
pairs=nchoosek(1:no_cells,2);
spike_times=cell(1,no_cells);

for n=1:no_cells
    
    spike_times{n}=t(spike_indices(cells==n));
    
end

for p=1:size(pairs,1)
    
    spike_times_1=spike_times{pairs(p,1)};
    spike_times_2=spike_times{pairs(p,2)};
    
    ISI_inter=[diff(spike_times_1) diff(spike_times_2)];
    
    if ~isempty(spike_times_1) && length(spike_times_1)<length(spike_times_2)
        
        spike_times_1((end+1):length(spike_times_2))=nan;
        
    elseif ~isempty(spike_times_2) && length(spike_times_2)<length(spike_times_1)
        
        spike_times_2((end+1):length(spike_times_1))=nan;
        
    end
    
    ISI_intra=nansum([spike_times_1; -spike_times_2]);
    
    ISIm(p)=mod(nanmean(ISI_intra)/nanmean(ISI_inter),1);
    
end



