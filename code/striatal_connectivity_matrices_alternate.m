function [CS,CG] = striatal_connectivity_matrices(no_cells);

	pct_both = .33;
	pct_gj_only = .18;
	pct_syn_only = .9;
	
	pct_gj = pct_both + pct_gj_only;
	pct_syn = pct_both + pct_syn_only;
	
	pct_gj_bi = .23;
	pct_syn_bi = .31;
	
	pairs = nchoosek(no_cells,2);
	
	gj_randos = rand(size(pairs,1),1);
	gj_indicator = randos < pct_gj;
	ngj_indicator = ~gj_indicator;
	
	syn_randos = rand(size(pairs,1),1);
	syn_only_indicator = ~gj_indicator & syn_randos < pct_syn_only/(1-pct_gj);
	both_indicator = gj_indicator & syn_randos < pct_both/pct_gj;
	syn_indicator = syn_only_indicator | both_indicator;
	
	gj_direction_randos = rand(size(pairs,1),1);
	gj_bi_indicator = gj_indicator & gj_direction_randos < pct_gj_bi;
	
	syn_direction_randos = rand(size(pairs,1),1);
	syn_bi_indicator = syn_indicator & syn_direction_randos < pct_syn_bi;
	
	reverse_pairs = fliplr(pairs);
	
	gj_pairs = [pairs(gj_indicator); reverse_pairs(gj_bi_indicator)];
	syn_pairs = [pairs(syn_indicator); reverse_pairs(syn_bi_indicator)];
	
end
	
	
	
	