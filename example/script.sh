epiphy-sim CpGs.bed sim.params -L sim_leaf.out -o sim_all.out
epiphy-est 7sp.nwk sim_leaf.out -v -o est_leaf.params
epiphy-est 7sp.nwk -c sim_all.out -v -o est_all.params
epiphy-post est_leaf.params sim_leaf.out -o est_leaf.post -v
epiphy-seg est_leaf.params est_leaf.post -o est_leaf.seg
