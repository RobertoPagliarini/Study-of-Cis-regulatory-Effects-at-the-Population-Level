function [AllChAEI,AEI_plusName,AEI_values,AllChAEI_Ent,MeanObsCisReg,MeanExpCisReg,MeanExpAEI,KLD_AEI,Entropy_AEI] = mainTestingAEI(Tissue)
%Main function of the tool to study Cis-regulatory Effects at the Population Level

AEI_plusName = {};

AEI_values = [];
 
AllChAEI = {};

AllChAEI_Ent = [];

ExpectedAEI_All = [];

N_ExpAEI_All = [];

ObservedCisRegulatoryDiversity_All = [];

N_ObsCisReg_All = [];

ExpectedCisRegulatoryDiversity_All = [];

N_ExpCisReg_All = [];


for i = 1:1:19

    disp(['Chromosome', num2str(i)]);

    [AEI,AEI_Ent,ExpectedAEI_Ch,N_ExpAEI_Ch,ObservedCisRegulatoryDiversity_Ch,N_ObsCisReg_Ch,ExpectedCisRegulatoryDiversity_Ch,N_ExpCisReg_Ch] = TestingAEI_v2(Tissue,i);

    AllChAEI{i,1} = AEI;

    AEI_values = [AEI_values;cell2mat(AllChAEI{i, 1}(2:end,8))];

    AllChAEI_Ent(i,1) = AEI_Ent; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 2:length(AEI)

        AEI_plusName{end+1,1} = AEI{j,1};

        AEI_plusName{end,2} = AEI{j,8};

        AEI_plusName{end,3} = AEI{j,9};

        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ExpectedAEI_All(i,1)  = ExpectedAEI_Ch;

    N_ExpAEI_All(i,1) = N_ExpAEI_Ch;

    ObservedCisRegulatoryDiversity_All(i,1) = ObservedCisRegulatoryDiversity_Ch;

    N_ObsCisReg_All(i,1) = N_ObsCisReg_Ch;

    ExpectedCisRegulatoryDiversity_All(i,1) = ExpectedCisRegulatoryDiversity_Ch;

    N_ExpCisReg_All(i,1) = N_ExpCisReg_Ch;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

MeanObsCisReg = sum(ObservedCisRegulatoryDiversity_All)/sum(N_ObsCisReg_All);

MeanExpCisReg = sum(ExpectedCisRegulatoryDiversity_All)/sum(N_ExpCisReg_All);

MeanExpAEI = sum(ExpectedAEI_All)/sum(N_ExpAEI_All);

KLD_AEI = getKullbackLeibler(ones(length(AEI_values),1)./length(AEI_values),AEI_values./sum(AEI_values));

Entropy_AEI = Entropy_Array(cell2mat(AEI(2:end,8)));

save(strcat('AEI_AllCh',Tissue,'.mat'), 'AllChAEI','AEI_plusName','AEI_values','AllChAEI_Ent','ExpectedAEI_All','N_ExpAEI_All','ObservedCisRegulatoryDiversity_All',...
     'N_ObsCisReg_All','ExpectedCisRegulatoryDiversity_All','N_ExpCisReg_All','MeanObsCisReg','MeanExpCisReg','MeanExpAEI','KLD_AEI','Entropy_AEI');