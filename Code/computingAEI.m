function [AEI,p,CisRegulatoryDiversity,N_CisReg] = computingAEI(AlleleExps,AllelesFrequencies)
%Computing allele expression imbalance from allele frequencies

HomozygousGenotypes = sum(AllelesFrequencies.^2);

HeterozygousGenotypes = 0;

CisRegulatoryDiversity = 0;

N_CisReg = length(AlleleExps);

for i = 1:length(AlleleExps)

    for j = i+1:length(AlleleExps)

        all_ij = [AlleleExps(i),AlleleExps(j)]; 

        n_ij = max(all_ij)/min(all_ij);
        
        HeterozygousGenotypes = HeterozygousGenotypes + 2*AllelesFrequencies(i)*AllelesFrequencies(j)*n_ij;

        CisRegulatoryDiversity = CisRegulatoryDiversity + log2(n_ij);

        N_CisReg = N_CisReg + 1;

    end

end

AEI = HomozygousGenotypes+HeterozygousGenotypes;

p = 1-HomozygousGenotypes;