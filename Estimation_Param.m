
%% CHEME 5999 Pset 3 - Parameter Estimations using Experimental Designs w/ Identifiable Parameters 
% The main script that estimates parameters
clear all
close all

EXP_NM = 3; % 3 out of 4 parameters identifiable with Experiment 3. No experiment for 4th parameter (Km1). 

t_i = 0;
t_f = 6000;   
t_inc = 300;  
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

randn('state',sum(100*clock));
rand('state',sum(100*clock));

while EXP_NM<=3
    %{
    if EXP_NM==4   
        [DF] = Opt_DataFile(DF,pset,P_Estimated_mtr,EXP_NM);
        Exp_Data_DF = Exp_Data(EXP_NM);
    %}
    if EXP_NM==3
        [DF] = DataFile (EXP_NM);
        Exp_Data_DF = Exp_Data(EXP_NM);
        Data = Exp_Data_DF.Data;
        Exp_max = max(Data);
        Exp_min = min(Data);
    elseif EXP_NM== 2
        [DF] = DataFile (EXP_NM);
        Exp_Data_DF = Exp_Data(EXP_NM);
        Data = Exp_Data_DF.Data;
        Exp_max = max(Data);
        Exp_min = min(Data);
    elseif EXP_NM==1
        [DF] = DataFile (EXP_NM);
        Exp_Data_DF = Exp_Data(EXP_NM);
        Data = Exp_Data_DF.Data;
        Exp_max = max(Data);
        Exp_min = min(Data);
        
    else disp ('Error,Construct level ranges from 1-5');
    end
    
    clearvars P_Estimated_mtr
    
    P = DF.Initial_Parameters;  

    [t,x] = Call_ODE(DF,tspan,P,EXP_NM); % The ODE function doesnt change with construct
    GP_t = transpose(x(:,4));

    for i = 1:(nstep+1)
        GP_t_Norm(i) = (GP_t(i) - GP_t(1))./(Exp_max - Exp_min); %normalize simulated data (GP_t) by experimental data (Exp_min/max) 
    end
    
    figure(EXP_NM)   %generate a profile with initial guess, pink line
    plot(t,GP_t_Norm,'m-','LineWidth',2); 
    xlabel('time'), ylabel('concentration')
    hold on

    GP_exp = Exp_Data_DF.Data;
    GP_exp_norm = GP_exp./(Exp_max);

    plot(t,GP_exp_norm,'b--'); %Experimental line
    hold on

    pset = DF.pset;
    
 %------------------------------------------------------
%This part of the program Estimates the identifiable parameters using the experimental
%data 

    loopsz = DF.iter_num;
    GP_opt_Norm_sum = zeros(1,(nstep+1));
    
    for i =1:(loopsz)
        [P_Estimated,P_all,GP_opt,fval,exitflag] = Estimation_fmincon(Exp_Data_DF,DF,tspan,pset,i);
        P_Estimated_mtr(i,:) = P_Estimated;
        
        if EXP_NM==3
            P_solution(i,:) = P_all; %Gathering final solution
        end
        
        for j=1:(nstep+1) %normalize GM_opt
            GP_opt_Norm(j) = (GP_opt(j)- GP_opt(1))./(Exp_max - Exp_min);       
        end
        
        GP_Mtx(i,:) = GP_opt_Norm;  
        GP_opt_Norm_sum = GP_opt_Norm_sum + GP_opt_Norm;
    end
    
    for j = 1:(nstep+1)
        GP_std(j) = std(GP_Mtx(:,j));  % find standard deviation of the plot
    end 

    GP_MEAN = GP_opt_Norm_sum./(loopsz); 
    GP_UB = GP_MEAN+1.96*GP_std; %95% confidence interval
    GP_LB = GP_MEAN-1.96*GP_std;
    
    plot(t,GP_Mtx ,'g');%individual modeled lines
    hold on
    plot(t,GP_MEAN,'k','LineWidth',2);% mean modeled line
    hold on
    plot(t,GP_UB,'k--',t,GP_LB,'k--');% 95% confidence
    hold all
    
    % --------------------------------------------------------
% part end

    if EXP_NM == 3 % all except 1  parameter estimated with Experiment 3 
        dlmwrite('P_solution2.txt',P_solution,'delimiter',' '); 
    end

    EXP_NM = EXP_NM + 1; % move on to the next construct

 end