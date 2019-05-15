%% CHEME 5999 Pset 3 -> ODEs Construct 
% Code contains ODE's for system 

function [t,x] = Call_ODE(DF,tspan,P,EXP_NM)
% This function solves sytem 

x0 = DF.Initial_Conditions; 

    function  dxdt = MassBalanceEqns(t,x) 
        %t,x,Km1,Km2,kcat,n) 
        %UG = x(1) ; Polypeptide = x(2); Mn2+ = x(3); GP = x(4) 
          % Km1 = P(1); Km2 = P(2); kcat = P(3); n = P(4);  
    F1 = 20; F2 = 20; F3 = 20; F4 = 20; V  = 40; % volume and flow rates 
    
        if EXP_NM == 1  % Experiment 1 conditions 
            Eo = 1.88; % [enzyme] in uM 
            x1_o = 5/Eo;  x2_o = 10/Eo; x3_o = 10/Eo; x4_o = 0/Eo; % initial concentrations
            theta = ((Eo*x(3))^P(4))/(50*Eo*x(3) + (Eo*x(3))^P(4));
            dxdt(1,1) = ((x1_o-x(1))*F1)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(2,1) = ((x2_o-x(2))*F2)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(3,1) = ((x3_o-x(3))*F3)/(P(3)*V) ; 
            dxdt(4,1) = theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2))) - (F4/(P(3)*V))*x(4) ; 

        elseif EXP_NM == 2  % Experiment 2 conditions
            Eo = 1.88;
            x1_o = 10/Eo;  x2_o = 5/Eo; x3_o = 10/Eo; x4_o = 0/Eo; % initial concentrations 
            theta = ((Eo*x(3))^P(4))/(50*Eo*x(3) + (Eo*x(3))^P(4));
            dxdt(1,1) = ((x1_o-x(1))*F1)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(2,1) = ((x2_o-x(2))*F2)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(3,1) = ((x3_o-x(3))*F3)/(P(3)*V) ; 
            dxdt(4,1) = theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2))) - (F4/(P(3)*V))*x(4) ; 
            
        elseif EXP_NM == 3  % Experiment 3 conditions
            Eo = 1.88;
            x1_o = 10/Eo;  x2_o = 10/Eo; x3_o = 5/Eo; x4_o = 0/Eo; % initial concentrations 
            theta = ((Eo*x(3))^P(4))/(50*Eo*x(3) + (Eo*x(3))^P(4));
            dxdt(1,1) = ((x1_o-x(1))*F1)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(2,1) = ((x2_o-x(2))*F2)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(3,1) = ((x3_o-x(3))*F3)/(P(3)*V) ; 
            dxdt(4,1) = theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2))) - (F4/(P(3)*V))*x(4) ; 
    
        elseif EXP_NM == 4  % Experiment 4 conditions
            Eo = 5.00;
            x1_o = 10/Eo;  x2_o = 10/Eo; x3_o = 10/Eo; x4_o = 0/Eo; % initial concentrations 
            theta = ((Eo*x(3))^P(4))/(50*Eo*x(3) + (Eo*x(3))^P(4));
            dxdt(1,1) = ((x1_o-x(1))*F1)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(2,1) = ((x2_o-x(2))*F2)/(P(3)*V) - theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2)));
            dxdt(3,1) = ((x3_o-x(3))*F3)/(P(3)*V) ; 
            dxdt(4,1) = theta*(x(1)/(P(1) + x(1)))*(x(2)/(P(2) + x(2))) - (F4/(P(3)*V))*x(4) ;
 
        end 
    end 
    
   [t,x]=ode45(@MassBalanceEqns,tspan,x0);
end 

