function [pval,ExpectedAEI,N_AEI,ExpectedCisRegulatoryDiversity,N_ExpCis,AEIpVec] = PermsStatTestAEI(AlleleExps,AllelesFrequencies, AEI)

pval = 0;

ExpectedAEI = 0;

N_AEI = 0;

ExpectedCisRegulatoryDiversity = 0;

AEIpVec = [];

N_ExpCis = length(AlleleExps);

if length(AlleleExps) < 7

     V = perms(AlleleExps);

     [nV,mV] = size(V);

     for i = 1:nV

         r = randperm(length(AllelesFrequencies));

         [AEI_p,~,CysRegDiv,N_CisReg] = computingAEI(V(i,:)',AllelesFrequencies(r));


         if AEI_p > AEI

             pval = pval+1;

         end

        ExpectedAEI = ExpectedAEI + AEI_p;

        ExpectedCisRegulatoryDiversity = ExpectedCisRegulatoryDiversity+CysRegDiv;

        N_ExpCis = N_ExpCis + N_CisReg;

        N_AEI = N_AEI + 1;

        AEIpVec = [AEIpVec;AEI_p];

     end

     pval = pval/nV;

else

    nperm = 1000;

    for i = 1:nperm

        r = randperm(length(AlleleExps));

        r2 = randperm(length(AlleleExps));

        Ap = AlleleExps(r);

        Af = AllelesFrequencies(r2);

        [AEI_p,~,CysRegDiv,N_CisReg] = computingAEI(Ap,Af);

         if AEI_p > AEI

             pval = pval+1;

         end

         ExpectedAEI = ExpectedAEI + AEI_p;

         ExpectedCisRegulatoryDiversity = ExpectedCisRegulatoryDiversity+CysRegDiv;

         N_ExpCis = N_ExpCis + N_CisReg;

         N_AEI = N_AEI + 1;

         AEIpVec = [AEIpVec;AEI_p];
         
    end

    pval = pval/nperm;

end

end