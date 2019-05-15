%% CHEME 5999 Problem Set 3 --> Experiment 1 "Fake" Measurements 
% Generate "fake" experimental data for Experiment 1 (x1o =[UGo] = 5/Eo & x1 = [UG] = 3/Eo); 

%clear all
%close all

global  x1_o x2_o x3_o x4_o F1 F2 F3 F4 V Eo 

Eo = 1.88; %uM

% set intitial conditions for state variables 
x1_o = 5/Eo; %uM 
x2_o = 10/Eo; %uM 
x3_o = 10/Eo ; %uM 
x4_o  = 0/Eo; 

% initial guess for S.S. value of state variables 
x1 = 3/Eo ; x2 = 5/Eo; x3 = 5/Eo; x4 = 3/Eo; 

x0 = [x1 x2 x3 x4] ; 

% parameters
F1 = 20; % uL/hr
F2 = 20; % uL/hr 
F4 = 20; % uL/hr
F3 = 20; %uL/hr 
V = 40; %uL 

% initial values 
% Original Parameter 
n = 4; 
kcat = 3600*1.5; % /hr
Km1 = 1.02; %uM    
Km2 = 2*Km1; %uM

%set time span
tspan = 0:5:100; 

% solve system of odes for original parameter set 
[t_orig,x_orig] = ode45(@(t,x) dxdt(t,x,Km1,Km2,kcat,n),tspan,x0); 

%Plot variation in each overtime
figure 
plot(t_orig,x_orig(:,1),'r');
hold on
plot(t_orig,x_orig(:,2),'b');
hold on
plot(t_orig,x_orig(:,3),'g');
hold on
plot(t_orig,x_orig(:,4),'m');
hold on
title ('EXP #1: Concentration vs. Time') 
xlabel('Time ');
ylabel('Concentration'); 
legend('dUG/dt', 'dP/dt','dM/dt','dGP/dt'); 
 


%save steady state values for each species 
x1_orig_ss = x_orig(end,1); x2_orig_ss = x_orig(end,2); x3_orig_ss = x_orig(end,3); x4_orig_ss = x_orig(end,4);
A = [' Steady State value of UG: ',num2str(x1_orig_ss)]; C = [' Steady State value of Mn2+: ',num2str(x3_orig_ss)];
B = [' Steady State value of Polypeptide: ',num2str(x2_orig_ss)]; D = [' Steady State value of GP: ',num2str(x4_orig_ss)];

disp(A); disp(B);  disp(C); disp(D);  

%table of experimental data
x1_points = x_orig(:,1).'; x2_points = x_orig(:,2).'; 
x3_points = x_orig(:,3).'; x4_points = x_orig(:,4).'; 
Data = [t_orig.' ; x1_points; x2_points; x3_points; x4_points ];
Data_table = Data.'; 

% write data into txt.file 
dlmwrite('Data_EXP1.txt',Data_table,'delimiter',' ','precision',5); 
type('Data_EXP1.txt') ; 


function func = dxdt(t,x,Km1,Km2,kcat,n)
global   F1 F2 F3 F4 V x1_o x2_o x3_o Eo  
% s1 = x(1) ; s2 = x(2); a = x(3); product = x(4) 

theta = ((Eo*x(3))^n)/(50*Eo*x(3) + (Eo*x(3))^n);

func1= ((x1_o-x(1))*F1)/(kcat*V) - theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2)));
func2= ((x2_o-x(2))*F2)/(kcat*V) - theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2)));
func3 = ((x3_o-x(3))*F3)/(kcat*V) ; 
func4 = theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))) - (F4/(kcat*V))*x(4) ; 

func = [func1; func2; func3; func4];
end 
