%% CHEME 5999 Pset 3 -> Solve Parameter Matrix (P_matrix) 

function [ P_mtrx ] = P_matrix(P,x,nstep,EXP_NM)
F1 = 20; F2 = 20; F3 = 20; F4 = 20; V  = 40; % volume and flow rates 

if EXP_NM == 1 
    Eo = 1.88; % [enzyme] in uM 
    x1_o = 5/Eo; x2_o = 10/Eo;  x3_o = 10/Eo; x4_o = 0/Eo;  % initial concentrations 

elseif EXP_NM == 2  
    Eo = 1.88; 
    x1_o = 10/Eo;  x2_o = 5/Eo; x3_o = 10/Eo; x4_o = 0/Eo; % initial concentrations

elseif EXP_NM == 3  
    Eo = 1.88; 
    x1_o = 10/Eo;  x2_o = 10/Eo; x3_o = 5/Eo; x4_o = 0/Eo; % initial concentrations

elseif EXP_NM == 4  
    Eo = 5.00; 
    x1_o = 10/Eo;  x2_o = 10/Eo; x3_o = 10/Eo; x4_o = 0/Eo; % initial concentrations

%{
if EXP_NM == 1 
    Eo = 1.88; % [enzyme] in uM 
    x1_o = 10/Eo; x2_o = 5/Eo;  x3_o = 20/Eo; x4_o = 0/Eo;  % initial concentrations 

elseif EXP_NM == 2  
    Eo = 4.00; 
    x1_o = 5/Eo;  x2_o = 5/Eo; x3_o = 3/Eo; x4_o = 0/Eo; % initial concentrations

elseif EXP_NM == 3  
    Eo = 1.88; 
    x1_o = 5/Eo;  x2_o = 10/Eo; x3_o = 10/Eo; x4_o = 0/Eo; % initial concentrations
 %}

end 

P_mtrx = zeros(4,4,(nstep+1)); % initialize parameter  matrix 

P_mtrx(3,1,:) = 0; 
P_mtrx(3,2,:) = 0; 
P_mtrx(3,4,:) = 0; 

for i = 2:(nstep+1)
     P_mtrx(1,1,i) = -( ((Eo*x(i,3))^P(4))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) )*( x(i,2)/(P(2) + x(i,2)) )*( x(i,1)/((P(1) + x(i,1))^2) );
     P_mtrx(1,2,i) = -( ((Eo*x(i,3))^P(4))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) )*( x(i,1)/(P(1) + x(i,1)) )*( x(i,2)/((P(2) + x(i,2))^2) );
     P_mtrx(1,3,i) = -( (x1_o - x(i,1))*F1 )/( V*((P(3))^2) ); 
     P_mtrx(1,4,i) = - ( x(i,1)/(P(1)+x(i,1)) )*( x(i,2)/(P(2)+x(i,2)) )*( (((Eo*x(i,3))^P(4))*log(Eo*x(i,3)))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) - ((Eo*x(i,3))^P(4))*((50*Eo*x(i,3) + (Eo*x(i,3))^P(4))^(-2))*(((Eo*x(i,3))^P(4))*log(Eo*x(i,3))) );
     P_mtrx(2,1,i) = -( ((Eo*x(i,3))^P(4))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) )*( x(i,2)/(P(2) + x(i,2)) )*( x(i,1)/((P(1) + x(i,1))^2) );
     P_mtrx(2,2,i) = -( ((Eo*x(i,3))^P(4))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) )*( x(i,1)/(P(1) + x(i,1)) )*( x(i,2)/((P(2) + x(i,2))^2) );
     P_mtrx(2,3,i) = -( (x2_o - x(i,2))*F2 )/( V*((P(3))^2) ); 
     P_mtrx(2,4,i) = - ( x(i,1)/(P(1)+x(i,1)) )*( x(i,2)/(P(2)+x(i,2)) )*( (((Eo*x(i,3))^P(4))*log(Eo*x(i,3)))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) - ((Eo*x(i,3))^P(4))*((50*Eo*x(i,3) + (Eo*x(i,3))^P(4))^(-2))*(((Eo*x(i,3))^P(4))*log(Eo*x(i,3))) );
     P_mtrx(3,3,i) = -( (x3_o - x(i,3))*F3 )/( V*((P(3))^2) );
     P_mtrx(4,1,i) = ( ((Eo*x(i,3))^P(4))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) )*( x(i,2)/(P(2) + x(i,2)) )*( x(i,1)/((P(1) + x(i,1))^2) );
     P_mtrx(4,2,i) = ( ((Eo*x(i,3))^P(4))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) )*( x(i,1)/(P(1) + x(i,1)) )*( x(i,2)/((P(2) + x(i,2))^2) );
     P_mtrx(4,3,i) = ( (x4_o - x(i,4))*F4 )/( V*((P(3))^2) ); 
     P_mtrx(4,4,i) = ( x(i,1)/(P(1)+x(i,1)) )*( x(i,2)/(P(2)+x(i,2)) )*( (((Eo*x(i,3))^P(4))*log(Eo*x(i,3)))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) - ((Eo*x(i,3))^P(4))*((50*Eo*x(i,3) + (Eo*x(i,3))^P(4))^(-2))*(((Eo*x(i,3))^P(4))*log(Eo*x(i,3))) );
     
end 
end 