function [pval, p, chi2,df, E] = hweChiSquare(geno)
% hweChiSquare - Chi-square test for Hardy-Weinberg equilibrium (k alleles)
%
% Input:
%   geno : k×k symmetric (or triangular) matrix of genotype counts
%
% Output:
%   chi2 : Chi-square statistic
%   pval : p-value
%   df   : degrees of freedom
%   E    : matrix of expected genotype counts under HWE
%
% 
% If pval > 0.05 → cannot reject HWE (genotypes consistent with allele frequencies).
% If pval < 0.05 → significant deviation → possible genotyping error, selection, inbreeding, or population structure.

% Example:
%   geno = [40 25 15; 0 20 10; 0 0 10];
%   [pval, chi2, df, E] = hweChiSquare(geno);

% -------------------------------------------------------------------------
% Author: 
% -------------------------------------------------------------------------


% --- Step 1: allele frequencies
x = geno;

N=sum(x(:)); %Total Genotypes

al=sum(x)+sum(x,2)'; %Alleles array

p = al / (2*N);

k = length(p);

% --- Step 2: expected genotype counts under HWE
E = zeros(k);
for i = 1:k
    for j = i:k
        if i == j
            E(i,j) = N * p(i)^2;
        else
            E(i,j) = 2 * N * p(i) * p(j);
        end
    end
end

% --- Step 3: chi-square statistic
chi2 = sum(((geno(triu(true(k))) - E(triu(true(k)))).^2) ./ E(triu(true(k))));

% --- Step 6: degrees of freedom and p-value
df = (k*(k+1)/2) - k;
pval = 1 - chi2cdf(chi2, df);
end