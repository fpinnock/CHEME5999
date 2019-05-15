%% CHEME 5999 Pset 3 -> Solve Sensitivity Matrix Analytically 

function [ SSM ] = Analytical_SSM( DF,x,nstep,t_inc,Prev_ID)

P = DF.Initial_Parameters;   %Obtain the initial guess vector from struct
P_size = DF.Num_Parameters;            
EXP_NM = DF.Construct;

[ J_mtrx ] = Jacobian( P,x,nstep,EXP_NM); % get the jacobian matrix
[ P_mtrx ] = P_matrix (P,x,nstep,EXP_NM ); % get the pmatrix
z0=[0;0;0;0];  % Initial value for Sensitivity (z) used in Ode45 

SSM = zeros((nstep),(P_size)); 

for i=2:(nstep+1)
    print=['solving for SSM at point ',num2str(i),'/',num2str(nstep+1)];
    disp(print);
    
    for j=1:(P_size)

        tspan=0:(t_inc):((i-1)*t_inc);  % to save computational time, only calculate z to the time point of interest
        [t,z]=ode15s(@senfunction,tspan,z0,[],J_mtrx,P_mtrx,t_inc,j); 
 
        SSM(i,j)=z(i,1);     % store z into SSM, here only take z(4), which is z_GP,j(t_k)
                             %(i.e. focus on sensitivity of experimental observable, glycoprotein product, to changes in each parameter)    
   
    end
end 

[SSM]= Normalize_SSM(DF,SSM,P,x,nstep);  % Identifiablity was estimated using an normalized SSM

if EXP_NM ==1
     
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0; %skips previous identified parameter and mark according columns 0
     end
    
elseif EXP_NM ==2

     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end
     
elseif EXP_NM ==3
    
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end
     
elseif EXP_NM ==4
   
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end

elseif EXP_NM ==5
   
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end
end 

end 

