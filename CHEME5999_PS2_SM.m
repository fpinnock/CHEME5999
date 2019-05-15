%% CHEME 5999 Problem Set 2, Question 2 
% Simulated annealing using the SA matlab algorithm (simulannealbnd)


% define upper and lower bound for each parameter (Km1, Km2, kcat, n) 
l_Km1 =.02; u_Km1 = 20; 
l_Km2 =.02; u_Km2 = 20;
l_kcat = 360*1.5; u_kcat = 36000*1.5; 
l_n = 0; u_n = 10;

lb = [l_Km1 l_Km2 l_kcat l_n];
ub = [u_Km1 u_Km2 u_kcat u_n];

% initial guess of parameter values
Km1_0 = 2; Km2_0 = 2;  kcat_0 = 3000; n_0 = 2; 

p0 = [Km1_0 Km2_0 kcat_0 n_0];

fun = @objective_fun;

options = optimoptions(@simulannealbnd,'FunctionTolerance',1e-6); 

p_simulannealbnd = simulannealbnd(fun,p0,lb,ub,options); 
disp(p_simulannealbnd); 

function f_error = objective_fun(p)
    %set time span
    t_final = 24; M = 10 ; 
    t_step = t_final/M; 
    tspan = 0:t_step:t_final ;
    
    % set initial conditions for each species   
    Eo = 1.88; %uM
    x1 = 5/Eo ; x2 = 3/Eo; x3 = 10/Eo; x4 = 3/Eo; 
    x0 = [x1 x2 x3 x4];  
    
    % define original parameter values and generate "fake" measurements
    Km1 = 1.02;  Km2 = 2*Km1; 
    kcat = 3600*1.5; n = 4; 
    [t_measured,x_measured] = ode45(@(t,x) dxdt(t,x,kcat,Km1,Km2,n),tspan,x0); 
    
    % now generate simulated data using parameter guesses
    Km1_guess = p(1); Km2_guess = p(2); kcat_guess = p(3); n_guess = p(4) ; 
    [t_simulated,x_simulated] = ode45(@(t,x) dxdt(t,x,Km1_guess,Km2_guess,kcat_guess,n_guess),tspan,x0);
    
    for i = 1:4 
        for j = 1:size(tspan)
           sqr_diff(j) = (x_simulated (i,j) - x_measured(i,j))^2;
           
        end
        summation_time(i) = sum(sqr_diff);
    end
    summation_species = sum(summation_time);  
  
    f_error = summation_species  
end

function func = dxdt(t,x,Km1,Km2,kcat,n) 
Eo = 1.88; %uM
% glycan (UG) = x(1) ; target polypeptide (P) = x(2); cofactor (Mn2+) = x(3); product (GP) = x(4) 
% x1_o = initial [UG], x2_o = intial [P], x3_o = initial [Mn2+], x4_o = intial [GP]

% define flow rates and volume 
F1 = 20; F2 = 20;  F3 = 20; F4 = 20; V  = 40;

% define dimensionless initial species concentrations
% x1_o = initial [UG], x2_o = intial [P]
% x3_o = initial [Mn2+], x4_o = intial [GP]
Eo = 1.88; %uM
x1_o = 10/Eo;  x2_o = 5/Eo; x3_o = 20/Eo ;  x4_o  = 0/Eo;


% enzyme activity dependency on cofactor w/ hill fxn  
theta = ((Eo*x(3))^n)/(50*Eo*x(3) + (Eo*x(3))^n); 

% define odes 
func1= (x1_o*F1)/(kcat*V) -  (x(1)*F4)/(kcat*V) - theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))); % dUG/dt
func2= (x2_o*F2)/(kcat*V) -  (x(2)*F4)/(kcat*V) - theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))); %dP/dt
func3 = (x3_o*F3)/(kcat*V) - (x(3)*F4)/(kcat*V) ; %dMn2+/dt
func4 = theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))) - (x(4)*F4)/(kcat*V) ;%dGP/dt 

func = [func1; func2; func3; func4];
end 