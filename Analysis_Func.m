%% CHEME 5999 Pset 3 - Analysis Function
% Reads estimated parameters and passes them to Analysis_Main.m to recreate
% parameterization trajectries.
function [DF] = Analysis_func(EXP_NM)

iter_num = 50;
filename = ['P_solution2.txt'];
file = fopen(filename);
P_est = fscanf(file, '%f',[iter_num,4]); 
     
%{
    if EXP_NM == 1 || 2 || 4    
         P_est(:,1)=0;
         P_est(:,2)=0;
         P_est(:,3)=0; 
         P_est(:,4)=0;
    end 
%}

    if EXP_NM == 1
        Eo = 1.88; % [enzyme] in uM 
        IC = [3/Eo          %UG_0 
              5/Eo          %Polypeptide_0
              5/Eo         %Mn2+_0 
              3/Eo          %GP_0
             ]; 
    elseif EXP_NM == 2  
        Eo = 1.88; % [enzyme] in uM 
        IC = [5/Eo          %UG_0 
              3/Eo          %Polypeptide_0
              5/Eo         %Mn2+_0 
              3/Eo          %GP_0
             ];  
         
    elseif EXP_NM == 3  
        Eo = 1.88; % [enzyme] in uM 
        IC = [5/Eo          %UG_0 
              5/Eo          %Polypeptide_0
              3/Eo         %Mn2+_0 
              3/Eo          %GP_0
             ];  
         
     elseif EXP_NM == 4  
        Eo = 5.00; % [enzyme] in uM 
        IC = [5/Eo          %UG_0 
              5/Eo          %Polypeptide_0
              5/Eo         %Mn2+_0 
              3/Eo          %GP_0
             ];  
    end 

DF.Num_Parameters = length(P_est(1,:));
DF.Initial_Conditions = IC;
DF.Construct = EXP_NM;
DF.Parameter_library = P_est;
DF.iter_num = iter_num;
end

    
    
    