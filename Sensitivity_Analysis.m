%% CHEME 5999 Pset 3 - Sensitivty Analysis to Design Parameterization Experiments (Step 1) 
% This program uses the initial guess parameter set to compute sensitivity
% matrices for each proposed construct. Then use "Identifiability.m" function
% to compute the identifiable paramters of each proposed construct

clear all
close all 

EXP_NM = input('Input Experiment Number:'); % Exp_Nm = "Experiment Number"  
Prev_ID = input('previously identified :');%(Include ALL previsouly identified parameter indices (single integers), input the array in a [], separated by ;)
Prev_ID = sort(Prev_ID);

% Set time step, t_k 
t_i = 0;
t_f = 6000;   
t_inc = 30; %Using finer time step for sensitivity matrix calculation, though not experimentally sampling at this fine of an interval  
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

if EXP_NM == 1
    DF = DataFile(EXP_NM);
    Exp_Data_DF= Exp_Data( EXP_NM );
elseif EXP_NM == 2
    DF = DataFile(EXP_NM);
    Exp_Data_DF= Exp_Data( EXP_NM );
elseif EXP_NM == 3
    DF = DataFile(EXP_NM);
    Exp_Data_DF= Exp_Data( EXP_NM );
elseif EXP_NM == 4
    DF = DataFile(EXP_NM);
    Exp_Data_DF= Exp_Data( EXP_NM );
elseif EXP_NM == 5
    DF = DataFile(EXP_NM);
    Exp_Data_DF= Exp_Data( EXP_NM );
else disp ('Error, Experiment numbers range from 1-5') 
    
end

P = DF.Initial_Parameters ; %obtain initial guess vector from "DataFile"  

[t,x] = Call_ODE(DF,tspan,P,EXP_NM); % The ODE function doesnt change along with construct

GP_t= transpose(x(:,4)); 

[ SSM ] = Analytical_SSM(DF,x,nstep,t_inc,Prev_ID);  %this function uses the analytical method to find SSM

[List,pset] = Identifiability (SSM);

%--------SSM heat map
A=abs(SSM') 
y=1:1.0:4;
x=0:5:100;
figure
imagesc(x,y,A)
caxis([0 0.5]) 
xlabel('Time(min)'), ylabel('Parameters)')
colorbar
fig=figure(1);
filename=['SSM',num2str(EXP_NM)];
print (fig,filename,'-dpng');
%--------

disp('Identifiable Parameters are:');
disp(pset);


