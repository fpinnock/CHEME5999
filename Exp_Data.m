%% CHEME 5999 Pset 3 -> Experimental Data 
% A function loads experimental data from text file and passes to main
% script

function [ Exp_Data_DF] = Exp_Data(EXP_NM)

if EXP_NM == 1
    Data= dlmread('Data_EXP1.txt'); 
    [r,c]=size(Data);
    time = Data(:,1);
    GP_prod = Data(:,2:c);  

    stepsz=size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GP_prod(i,:));
        lowbound(i)=min(GP_prod(i,:));
        avg(i)=sum(GP_prod(i,:))/(c-1); 
    end
    
elseif EXP_NM == 2
    Data = dlmread('Data_EXP2.txt'); 
    [r,c] = size(Data);
    time = Data(:,1);
    GP_prod = Data(:,2:c);  

    stepsz = size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GP_prod(i,:));
        lowbound(i)=min(GP_prod(i,:));
        avg(i)=sum(GP_prod(i,:))/(c-1); 
    end 

elseif EXP_NM == 3
    Data = dlmread('Data_EXP3.txt'); 
    [r,c] = size(Data);
    time = Data(:,1);
    GP_prod = Data(:,2:c);  

    stepsz = size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GP_prod(i,:));
        lowbound(i)=min(GP_prod(i,:));
        avg(i)=sum(GP_prod(i,:))/(c-1); 
    end 
    
 elseif EXP_NM == 4
    Data = dlmread('Data_EXP4.txt'); 
    [r,c] = size(Data);
    time = Data(:,1);
    GP_prod = Data(:,2:c);  

    stepsz = size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GP_prod(i,:));
        lowbound(i)=min(GP_prod(i,:));
        avg(i)=sum(GP_prod(i,:))/(c-1); 
    end 
    
end 

Exp_Data_DF.Data = GP_prod(:,4);
Exp_Data_DF.avg = avg';
Exp_Data_DF.Upbound = Upbound';
Exp_Data_DF.Lowbound = lowbound';
Exp_Data_DF.timestep = stepsz;

end 

