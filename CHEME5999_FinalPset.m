% CHEME5999 Final Assignment 

clear all
close all 
%% Model Simulation 

global  x1_o x2_o x3_o x4_o F1 F2 F3 F4 V Eo 

Eo = 1.88; %uM

% initial conditions for dimensionless variables 
% (x1 = UG   x2= P  x3 = Mn2+  x4 = GP) 
x1_o = 10/Eo; x2_o = 5/Eo;  x3_o = 20/Eo ; x4_o  = 0/Eo; 

% initial guesses
x1 = 5/Eo ; x2 = 3/Eo; x3 = 10/Eo; x4 = 3/Eo; 
x0 = [x1 x2 x3 x4] ; 

% parameters
F1 = 20; F2 = 20;  F4 = 20; F3 = 20; %uL/hr 
V = 40; %uL 

% original parameter set  
n = 4; kcat = 3600*1.5; %/hr
Km1 = 1.02;  Km2 = 2*Km1; %uM


%set time span
tspan = [0 24]; 

% solve system of odes for model simulation 
[t_mod,x_mod] = ode45(@(t,x) dxdt(t,x,Km1,Km2,kcat,n),tspan,x0);


%% Model Simulation vs. Fake Experimental Data with "Real" Parameters 

% Real Parameters (Gerber, S., et al, J. Biol. Chem. 2013,288, 8849-8861)  
n_real = 8.8; kcat_real = 5912.2; % /hr
Km1_real = 4.58 ; Km2_real = 0.02; % or (1.02 +/- 0.06 uM , 7.88+/- 1.76 mM)   

% solve system of odes for  real parameter set  
[t_real,x_real] = ode45(@(t,x) dxdt(t,x,Km1_real,Km2_real,kcat_real,n_real),tspan,x0); 

% Plot Model Simulation vs. Real Parameters
figure 
plot(t_mod,x_mod(:,1),'r','Linewidth',2);
hold on
plot(t_mod,x_mod(:,2),'b','Linewidth',2);
hold on
plot(t_mod,x_mod(:,3),'g','Linewidth',2);
hold on
plot(t_mod,x_mod(:,4),'m','Linewidth',2);
hold on
plot(t_real,x_real(:,1),'--r','Linewidth',2);
hold on
plot(t_real,x_real(:,2),'--b','Linewidth',2);
hold on
plot(t_real,x_real(:,3),'--g','Linewidth',2);
hold on
plot(t_real,x_real(:,4),'--m','Linewidth',2);
title ('Model Simulation vs. Real Parameters ') 
xlabel('Time ');
ylabel('Concentration'); 
legend('dUG/dt Model', 'dP/dt Model','dM/dt Model','dGP/dt Model','dUG/dt real', 'dP/dt real','dM/dt real','dGP/dt real'); 


%% Model Simulation vs. Fake Experimental Data with Learned Parameters from Simulated Annealing Algorithm (SM) 

% SM parameters (see "CHEME5999_PS2_SM " for code) 
n_SM = 0; kcat_SM = 2804.8; %/hr
Km1_SM = 13.2;  Km2_SM = 19.9; %uM (13.6 18.7 2978.7 0) 

% solve system of odes for model simulation 
[t_SM,x_SM] = ode45(@(t,x) dxdt(t,x,Km1_SM,Km2_SM,kcat_SM,n_SM),tspan,x0);

% Plot Model Simulation vs. SM Experimental Data 
figure 
plot(t_mod,x_mod(:,1),'r','Linewidth',2);
hold on
plot(t_mod,x_mod(:,2),'b','Linewidth',2);
hold on
plot(t_mod,x_mod(:,3),'g','Linewidth',2);
hold on
plot(t_mod,x_mod(:,4),'m','Linewidth',2);
hold on
plot(t_SM,x_SM(:,1),'--r','Linewidth',2);
hold on
plot(t_SM,x_SM(:,2),'--b','Linewidth',2);
hold on
plot(t_SM,x_SM(:,3),'--g','Linewidth',2);
hold on
plot(t_SM,x_SM(:,4),'--m','Linewidth',2);
title ('Model Simulation vs. Simulated Annealing Parameters ') 
xlabel('Time ');
ylabel('Concentration'); 
legend('dUG/dt Model', 'dP/dt Model','dM/dt Model','dGP/dt Model','dUG/dt SM', 'dP/dt SM','dM/dt SM','dGP/dt SM'); 

%% Model Simulations vs. Fake Experimental Data w/ Learned Parameters using Experimental Design Approach 
% See folder entitled "CHEME5999_ED_Algorithm" 

% Experimental Deisgn(ED) Parameter Sets 
n_ED = 6.5354; kcat_ED = 168670; 
Km1_ED = 0; Km2_ED = 3.1459;

% solve system of odes for  simulated parameter set  
[t_ED,x_ED] = ode45(@(t,x) dxdt(t,x,Km1_ED,Km2_ED,kcat_ED,n_ED),tspan,x0);

% Plot Model Simulation vs.Experimental Design Parameters ( see "CHEME5999_PS3" folder)  
figure 
plot(t_mod,x_mod(:,1),'r','Linewidth',2);
hold on
plot(t_mod,x_mod(:,2),'b','Linewidth',2);
hold on
plot(t_mod,x_mod(:,3),'g','Linewidth',2);
hold on
plot(t_mod,x_mod(:,4),'m','Linewidth',2);
hold on
plot(t_ED,x_ED(:,1),'--r','Linewidth',2);
hold on
plot(t_ED,x_ED(:,2),'--b','Linewidth',2);
hold on
plot(t_ED,x_ED(:,3),'--g','Linewidth',2);
hold on
plot(t_ED,x_ED(:,4),'--m','Linewidth',2);
title ('Model Simulation vs. Experimental Design Parameters ') 
xlabel('Time ');
ylabel('Concentration'); 
legend('dUG/dt Model', 'dP/dt Model','dM/dt Model','dGP/dt Model','dUG/dt ED', 'dP/dt ED','dM/dt ED','dGP/dt ED'); 


%%  Functions 

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