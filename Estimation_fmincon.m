%% CHEME 5999 Pset 3 - Estimate Identifiable Parameters using fmincon
% A function that estimates identifiable parameters by minimizing errors to experimental data 

function [P_Estimated,P_all,GP_opt,fval,exitflag]=Estimation_fmincon(Exp_Data_DF,DF,tspan,pset,i)

    EXP_NM = DF.Construct;
    P_lib = DF.Parameter_library;
    NParameters = DF.Num_Parameters;
    iter_num = DF.iter_num;
    P = P_lib(i,:); %Columns are the different parameters, rows are different sets (1000 sets total) 

    function [Error] = fminconfunc(y)

       timesize = Exp_Data_DF.timestep; % get experimental data 

       for k = 1:length(pset)  %update P with identifiable parameters
          P(pset(k)) = y(k); % define identifiable parameter 
       end

       [t,x] = Call_ODE(DF,tspan,P,EXP_NM); % solve ODE with original guess for identifiable parameter in designated EXP 
       GP_est = transpose(x(:,4)); %Basing estimate error on GFP prediction
       Error = 0; 

       for j=1:(timesize)
          Error = Error + (GP_est(j)- Exp_Data_DF.avg(j)).^2; %Sum of Error 
       end

    end

   y0 = P(pset(1));
    
   [A,b] = fmincon_constraints(pset,P);
   for m = 2:length(pset)
        y0 = [y0,P(pset(m))];
   end

   lb = 0.33*y0; 
   ub = 3.33*y0;
   
   [y,fval,exitflag] = fmincon(@fminconfunc,y0,A,b,[],[],lb,ub);
   P_Estimated = y;

    P_all = P;
    for i = m:(length(pset))
        P_all((pset(m)))=P_Estimated(m);
    end
    
    P = P_all;
    [t,x] = Call_ODE(DF,tspan,P,EXP_NM);
    GP_opt=transpose(x(:,4));    
    
end
   