%% CHEME 5999 Pset 3 -> Solve Jacobian Matrix 

function [ J_mtrx ] = Jacobian(P,x,nstep,EXP_NM)

F1 = 20; F2 = 20; F3 = 20; F4 = 20; V  = 40; % volume and flow rates

if EXP_NM == 1
    Eo = 1.88; % [enzyme] in uM 
elseif EXP_NM == 2 
    Eo = 1.88;
elseif EXP_NM == 3 
    Eo = 1.88;
elseif EXP_NM == 4 
    Eo = 5.00;
end 

J_mtrx = zeros(4,4,(nstep+1)); % initialize Jacobian matrix 

J_mtrx(1,4,:) = 0;  % These elements are constant 
J_mtrx(2,4,:) = 0; 
J_mtrx(3,1,:) = 0;  
J_mtrx(3,2,:) = 0; 
J_mtrx(3,3,:) = -F3/(P(3)*V) ;
J_mtrx(3,4,:) = 0; 
J_mtrx(4,4,:) = -F4/(P(3)*V) ; 


for i = 2:(nstep+1)
   J_mtrx(1,1,i) = -F1/(P(3)*V) -  ( (Eo*x(i,3))^P(4) )/( 50*Eo*x(i,3) + (Eo*x(i,3))^P(4) )*( x(i,2)/(P(2)+x(i,2)) )*( (P(1)+x(i,1))^(-1) - x(i,1)*(P(1) + x(i,1))^(-2) ) ;
   J_mtrx(1,2,i) = -( (Eo*x(i,3))^P(4) )/( 50*Eo*x(i,3) + (Eo*x(i,3))^P(4) )*( x(i,1)/(P(1)+x(i,1)) )*( (P(2)+x(i,2))^(-1) - x(i,2)*(P(2) + x(i,2))^(-2) ) ;
   J_mtrx(1,3,i) = -( x(i,1)/(P(1)+x(i,1)) )*( x(i,2)/(P(2)+x(i,2)) )*( (P(4)*Eo*(Eo*x(i,3)^(P(4)-1))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) ) -  ((50*Eo*x(i,3) + (Eo*x(i,3))^P(4))^(-2))*(50*Eo + P(4)*Eo*(Eo*x(i,3))^(P(4)-1))*(Eo*x(i,3))^P(4) );
   J_mtrx(2,1,i) = -( (Eo*x(i,3))^P(4) )/( 50*Eo*x(i,3) + (Eo*x(i,3))^P(4) )*( x(i,2)/(P(2)+x(i,2)) )*( (P(1)+x(i,1))^(-1) - x(i,1)*(P(1) + x(i,1))^(-2) ) ;
   J_mtrx(2,2,i) = - F2/(P(3)*V) -  ( (Eo*x(i,3))^P(4) )/( 50*Eo*x(i,3) + (Eo*x(i,3))^P(4) )*( x(i,1)/(P(1)+x(i,1)) )*( (P(2)+x(i,2))^(-1) - x(i,2)*(P(2) + x(i,2))^(-2) ) ;
   J_mtrx(2,3,i) = -( x(i,1)/(P(1)+x(i,1)) )*( x(i,2)/(P(2)+x(i,2)) )*( (P(4)*Eo*(Eo*x(i,3)^(P(4)-1))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) ) -  ((50*Eo*x(i,3) + (Eo*x(i,3))^P(4))^(-2))*(50*Eo + P(4)*Eo*(Eo*x(i,3))^(P(4)-1))*(Eo*x(i,3))^P(4) );
   J_mtrx(4,1,i) = ( (Eo*x(i,3))^P(4) )/( 50*Eo*x(i,3) + (Eo*x(i,3))^P(4) )*( x(i,2)/(P(2)+x(i,2)) )*( (P(1)+x(i,1))^(-1) - x(i,1)*(P(1) + x(i,1))^(-2) ) ;
   J_mtrx(4,2,i) = ( (Eo*x(i,3))^P(4) )/( 50*Eo*x(i,3) + (Eo*x(i,3))^P(4) )*( x(i,1)/(P(1)+x(i,1)) )*( (P(2)+x(i,2))^(-1) - x(i,2)*(P(2) + x(i,2))^(-2) ) ; 
   J_mtrx(4,3,i) = ( x(i,1)/(P(1)+x(i,1)) )*( x(i,2)/(P(2)+x(i,2)) )*( (P(4)*Eo*(Eo*x(i,3)^(P(4)-1))/(50*Eo*x(i,3) + (Eo*x(i,3))^P(4)) ) -  ((50*Eo*x(i,3) + (Eo*x(i,3))^P(4))^(-2))*(50*Eo + P(4)*Eo*(Eo*x(i,3))^(P(4)-1))*(Eo*x(i,3))^P(4) );
  
end 