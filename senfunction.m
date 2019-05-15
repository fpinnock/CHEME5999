%% CHEME 5999 Pset 3 -> Solve for Sensitivity Matrix 

% Sensitivity coefficient matrix , z_ij(t), is given by a set of ODEs
% see Methods section of (DOI: 10.1021/acssynbio.5b00077) 

function [DZDT] = senfunction(t,z,J_mtrx,P_mtrx,t_inc,j) %M,PM)
 
        k = round(1+t./t_inc) ;
        JM = J_mtrx(:,:,k); % every time pass a 4x4 jacbian matrix to the ODE function
        PM = P_mtrx(:,j,k); % every time pass a 4x1 pmatrix to the ODE function   
        
        DZDT = JM*z + PM;

end 