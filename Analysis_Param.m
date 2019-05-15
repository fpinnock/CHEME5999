%% CHEME 5999 Pset 3 - Analysis 
% Recreates parameterization experiment trajectories using estimated
% parameters
% ------------------------------------------------------------------------------------- %

clear all
close all
EXP_NM = 3; % again, only care about Experiment 3 output  

t_i = 0;
t_f = 6000;   
t_inc = 300;  
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

while EXP_NM <= 3
    %{
    if EXP_NM == 4    
        [DF] = Analysis_Func(EXP_NM);
        Exp_Data_DF= Exp_Data(EXP_NM);
    %}
     if EXP_NM == 3
        [DF] = Analysis_Func(EXP_NM);
        Exp_Data_DF = Exp_Data( EXP_NM );
        Data = Exp_Data_DF.Data;
        Exp_max = max(Data);
        Exp_min = min(Data);
    elseif EXP_NM == 2
        [DF] = Analysis_Func(EXP_NM);
        Exp_Data_DF = Exp_Data( EXP_NM );
        Data = Exp_Data_DF.Data;
        Exp_max = max(Data);
        Exp_min = min(Data);
    elseif EXP_NM == 1
        [DF] = Analysis_Func(EXP_NM);
        Exp_Data_DF = Exp_Data( EXP_NM );
        Data = Exp_Data_DF.Data;
        Exp_max = max(Data);
        Exp_min = min(Data);
    end

    Norm_Exp_Avg = Exp_Data_DF.avg./Exp_max;
    iter_num = DF.iter_num;
    loopsz = iter_num;
    GP_opt_Norm_sum = zeros(1,(nstep+1));
    GP_error = zeros(iter_num,1);
    
    for i = 1:(loopsz)
        P_lib = DF.Parameter_library; % P_lib = P_est = P_solution2.txt
        NParameters = DF.Num_Parameters;
        
        %{
        for j = 1:(NParameters)        
          index = randi(loopsz,1); % randomly choose estimated parameters from P_solution2.txt
          P(j) = P_lib(index,j);
        end
        %}
        index = randi(loopsz,1); 
        P = P_lib(index,:); % use estimated parameter set from randomly selected row of P_solutio2.txt
        
        [t,x] = Call_ODE(DF,tspan,P,EXP_NM);
        GP_opt = transpose(x(:,4)); 
        for j = 1:(nstep+1)
            GP_opt_Norm(j) = (GP_opt(j) - GP_opt(1))./(Exp_max - Exp_min);  
        end
    
        GP_Mtx(i,:) = GP_opt_Norm ; 
        GP_opt_Norm_sum = GP_opt_Norm_sum + GP_opt_Norm;
            
    end
    
    for j=1:(nstep+1)
        GP_std(j)=std(GP_Mtx(:,j));
    end 


    GP_MEAN = GP_opt_Norm_sum./(loopsz);
    GP_UB = GP_MEAN + 1.96*GP_std;
    GP_LB = GP_MEAN - 1.96*GP_std;
    GP_exp = Exp_Data_DF.Data;
    GP_exp_norm = GP_exp./(Exp_max);

    figure(EXP_NM)   %generate a profile 
    %plot(t,GP_Mtx ,'g');
    %axis([0 6000 0 1]);
    %hold all
    plot(t,GP_MEAN,'k','LineWidth',3);
    hold on 
    plot(t,GP_UB,'k--')
    hold on
    plot(t,GP_LB,'k--');
    hold on 
    plot(t,GP_exp_norm,'b--');
    xlabel('time'); 
    ylabel('concentration');  
    legend('GP Normalized Optimal', 'GP Mean','GP Lower Bound','GP Upper Bound', 'GP Normalized Exerpimental '); 
    
    hold on 
    DATA=[GP_MEAN;GP_UB;GP_LB;GP_exp_norm']; % Saved text files contain mean,lower bound, upper bound and normalized experimenl data.
    filename=['DATA_EXP_',num2str(EXP_NM),'.txt'];
    dlmwrite(filename,DATA','delimiter',' ');
   
    EXP_NM = EXP_NM + 1;
end

