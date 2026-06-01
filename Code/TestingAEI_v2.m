function [AEI,AEI_Ent,ExpectedAEI_Ch,N_ExpAEI_Ch,ObservedCisRegulatoryDiversity_Ch,N_ObsCisReg_Ch,ExpectedCisRegulatoryDiversity_Ch,N_ExpCisReg_Ch] = TestingAEI_v2(Tissue,Chromosome)

%Selecting chromosome
addpath(strcat('/Users//Library/CloudStorage/OneDrive-UniversitàdegliStudidiUdine/Ricerca scientifica/Serie Alleliche/Data/all_genes_nofilter-2/chr',num2str(Chromosome),'/output'));

%Selecting tissue
if strcmp(Tissue,'Leaves')

    addpath(strcat('/Users//Ricerca scientifica/Serie Alleliche/Allele Imbalance/LeavesResultsDataMat/Chromosome',num2str(Chromosome),'AlleleValues/'))

    Allres = dir(strcat('/Users//Ricerca scientifica/Serie Alleliche/Allele Imbalance/LeavesResultsDataMat/Chromosome',num2str(Chromosome),'AlleleValues/*.mat'));

end

if strcmp(Tissue,'BerriesSoft')

    addpath(strcat('/Users//Ricerca scientifica/Serie Alleliche/Allele Imbalance/BerriesSoftResultsDataMat/Chromosome',num2str(Chromosome),'AlleleValues/'))

    Allres = dir(strcat('/Users//Ricerca scientifica/Serie Alleliche/Allele Imbalance/BerriesSoftResultsDataMat/Chromosome',num2str(Chromosome),'AlleleValues/*.mat'));

end

if strcmp(Tissue,'BerriesHard')

    addpath(strcat('/Users//Ricerca scientifica/Serie Alleliche/Allele Imbalance/BerriesHardResultsDataMat/Chromosome',num2str(Chromosome),'AlleleValues/'))

    Allres = dir(strcat('/Users//Ricerca scientifica/Serie Alleliche/Allele Imbalance/BerriesHardResultsDataMat/Chromosome',num2str(Chromosome),'AlleleValues/*.mat'));

end

%Setting output
AEI{1,1} = 'GeneName';

AEI{1,2} = 'AlleleNames';

AEI{1,3} = 'AlleleResults';

AEI{1,4} = 'AlleleFrequencies';

AEI{1,5} = 'FrequencieErrors';

AEI{1,6} = 'ObsGenotypes';

AEI{1,7} = 'GenotypesFrequencies';

AEI{1,8} = 'AEIAlleles';

AEI{1,9} = 'AEIAlleles_permpval';

AEI{1,10} = 'AEIAlleles_heteroprobability';

AEI{1,11} = 'AEIpermutation';

AEI{1,12} = 'TestingHWE';

AEI{1,13} = 'TotalGenotypes';

AEI{1,14} = 'HWEChiSquare';

AEI{1,15} = 'HWEChiSquareEGC';

ObservedCisRegulatoryDiversity_Ch = 0;

N_ObsCisReg_Ch = 0;

ExpectedCisRegulatoryDiversity_Ch = 0;

N_ExpCisReg_Ch = 0;

ExpectedAEI_Ch = 0;

N_ExpAEI_Ch = 0;

%AEI computing
for i = 1:length(Allres)

   Allele_Results = [];

   AlleleNames = {};

   %File name
   FileName = Allres(i).name(1:end-4);

   data = load(Allres(i).name);

   data = data.prova;

   [m,n] = size(data);

   %Saving allele names and allele results
   for ii = 2:m

       Allele_Results(ii-1,1) = data{ii,2};

       AlleleNames{ii-1,1} = data{ii,1};

   end


   %Observed genotypes and allele frequencies
   [ObsGenotype,fa,e_fa,fg,al,Ng] = readObsGenotype(FileName,AlleleNames);

    if isempty(ObsGenotype) == 0

        AEI{i+1,1} = FileName;

        AEI{i+1,2} =  AlleleNames;

        AEI{i+1,3} = Allele_Results;

        AEI{i+1,4} = fa';

        AEI{i+1,5} = e_fa;

        AEI{i+1,6} = ObsGenotype;

        AEI{i+1,7} = fg;

        %Computing allele expression imbalance from allele frequencies
        [aei,p,ObservedCisRegulatoryDiversity_g,N_ObsCisReg_g] =  computingAEI(Allele_Results,fa');

        AEI{i+1,8} =  aei;

        [pval,ExpectedAEI_g,N_ExpAEI_g,ExpectedCisRegulatoryDiversity_g,N_ExpCisReg_g,AEIpVec] = PermsStatTestAEI(Allele_Results,fa', aei);
        
        AEI{i+1,9} =  pval;

        AEI{i+1,10} = p;

        %Computing allele expression imbalance from genotypes 
        AEI{i+1,11} = AEIpVec;
        
        %Testing Hardy-Weinberg
        AEI{i+1,12} = hwetest(ObsGenotype);

        %Total genotypes
        AEI{i+1,13} = Ng;

        ObservedCisRegulatoryDiversity_Ch = ObservedCisRegulatoryDiversity_Ch + ObservedCisRegulatoryDiversity_g;

        N_ObsCisReg_Ch = N_ObsCisReg_Ch + N_ObsCisReg_g ;

        ExpectedCisRegulatoryDiversity_Ch = ExpectedCisRegulatoryDiversity_Ch + ExpectedCisRegulatoryDiversity_g;

        N_ExpCisReg_Ch = N_ExpCisReg_Ch + N_ExpCisReg_g;

        ExpectedAEI_Ch = ExpectedAEI_Ch + ExpectedAEI_g;

        N_ExpAEI_Ch = N_ExpAEI_Ch + N_ExpAEI_g;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AEI{i+1,14} = hweChiSquare(ObsGenotype);

        AEI{i+1,15} = hweChiSquareEGC(ObsGenotype,Allele_Results);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    end

end

[AEI_Ent] = Entropy_Array(cell2mat(AEI(2:end,8)));



