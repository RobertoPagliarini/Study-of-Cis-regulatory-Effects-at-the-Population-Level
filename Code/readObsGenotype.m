function [ObsGenotype,fa,e_fa,fg,al,Ng] = readObsGenotype(FileName, WName)

ObsGenotype = zeros(length(WName),length(WName));

A = readtable(strcat(FileName,'_serie_allelica.txt'));

[n,m] = size(A);

for i = 1:n

    index_hap = [];

    for j = 1:length(WName)

        if strcmp(A.hap_A{i,1}, WName{j,1})

            index_hap(end+1) = j;

        end

        if strcmp(A.hap_B{i,1}, WName{j,1})

            index_hap(end+1) = j;

        end

    end

    if length(index_hap) == 2

        s = sort(index_hap);

        ObsGenotype(s(1),s(2)) = ObsGenotype(s(1),s(2)) + 1;

    end

end

%Total Genotypes
Ng=sum(ObsGenotype(:));

%Alleles array
al=sum(ObsGenotype)+sum(ObsGenotype,2)'; 

%Allelic frequencies
fa=al./(2*Ng);

%Error of frequencies determination
e_fa=sqrt(prod(fa)/Ng); 

mfa=(fa'*fa); %matrix of binomial expansion of allelic frequencies

xe=round(mfa.*Ng); %Expected matrix under HWP

%Matrix genotypes frequencies
fg=ObsGenotype./Ng;
