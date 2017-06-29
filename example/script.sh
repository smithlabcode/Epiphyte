##### Small test examples #####
# Files coming with the package: CpGs.bed  sim.params 7sp.nwk
# Modify the following line to go to the "example" subdirectory of Epiphyte
cd <path to Epiphyte>/example

# simulate methylomes 
../bin/epiphy-sim CpGs.bed sim.params -L sim_leaf.out -o sim_all.out
# estimate model parameters from leaf species observations
../bin/epiphy-est 7sp.nwk sim_leaf.out -v -o est_leaf.params
# estimate model parameters with complete observations
../bin/epiphy-est 7sp.nwk -c sim_all.out -v -o est_all.params
# get posterior methylation probabilities 
../bin/epiphy-post est_leaf.params sim_leaf.out -o est_leaf.post -v
# get methylation state segmentations
../bin/epiphy-seg est_leaf.params est_leaf.post -o est_leaf.seg
