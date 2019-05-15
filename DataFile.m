%% CHEME 5999 Pset 3 -> Data File 

function [DF] = DataFile (EXP_NM) 

iter_num = 50; 
pset = 3; % identifiable parameters for EXP_NM = #1 
% Set initial guess for S.S. values of state variables (species concentrations)  
if EXP_NM == 1
    Eo = 1.88; % [enzyme] in uM 
    IC = [3/Eo          %UG 
          5/Eo          %Polypeptide
          5/Eo         %Mn2+ 
          3/Eo          %GP
         ]; 
elseif EXP_NM == 2
    Eo = 1.88; % [enzyme] in uM 
    IC = [5/Eo          %UG 
          3/Eo          %Polypeptide
          5/Eo         %Mn2+
          3/Eo          %GP
         ];  
elseif EXP_NM == 3  
    Eo = 1.88; % [enzyme] in uM 
    IC = [5/Eo          %UG
          5/Eo          %Polypeptide
          3/Eo         %Mn2+ 
          3/Eo          %GP
         ];  
elseif EXP_NM == 4  
    Eo = 5.00; % [enzyme] in uM 
    IC = [5/Eo          %UG 
          5/Eo          %Polypeptide
          5/Eo         %Mn2+ 
          3/Eo          %GP
         ];  
 
end 

%{
if EXP_NM == 1
    Eo = 1.88; % [enzyme] in uM 
    IC = [5/Eo          %UG_0 
          3/Eo          %Polypeptide_0
          10/Eo         %Mn2+_0 
          3/Eo          %GP_0
         ]; 
elseif EXP_NM == 2
    Eo = 4.00; % [enzyme] in uM 
    IC = [3/Eo          %UG_0 
          3/Eo          %Polypeptide_0
          1/Eo         %Mn2+_0 
          3/Eo          %GP_0
         ];  
elseif EXP_NM == 3  
    Eo = 1.88; % [enzyme] in uM 
    IC = [3/Eo          %UG_0 
          5/Eo          %Polypeptide_0
          7/Eo         %Mn2+_0 
          3/Eo          %GP_0
         ];  
end 
%} 


%Initial guess for all parameters (Km1, Km2, kcat,n); guess the same for
%both Experiments and based on PSET 1 Literature Values 
P_i = [ 0    %KM1 
        3   %KM2 
        54000  %kcat
        7];       %n 

P_i_lb = P_i.*0.85;    %artificial range, within 15% of the best guess
P_i_ub = P_i.*1.15; 


NParameters = length(P_i);
Nstates = length(IC);


a = P_i_lb;      %generate a set of P from 15% range about primary guess, first round, all random
 b = P_i_ub ;
 for j=1:(iter_num)
     for i=1:(NParameters)
         P(j,i)=a(i)+(b(i)-a(i)).*rand(1,1);
     end
 end

DF.Initial_Parameters = P_i;
DF.Initial_Conditions = IC;
DF.Num_States = Nstates;
DF.Num_Parameters = NParameters;
DF.Construct = EXP_NM;
DF.P_i_LB = P_i_lb;
DF.P_i_UB = P_i_ub;
DF.Parameter_library = P;
DF.iter_num = iter_num;
DF.pset=pset;
end 

