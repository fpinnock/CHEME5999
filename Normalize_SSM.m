%% CHEME 5999 Pset 3 -> Normalize SSM

function [SSM] = Normalize_SSM (DF,SSM,P,x,nstep)
c = 2;
NSSM = zeros(nstep./10,(DF.Num_Parameters));

for i=11:10:(nstep+1) % Pulling out coarse-grained time steps that corresponds to experiments
    
    for j=1:(DF.Num_Parameters)
        NSSM(c,j)= SSM(i,j).*(P(j)./x(i,1)); % remember, only focused on x(4) corresponding to glycoprotein product GP
    end
    c=c+1;
end

SSM=NSSM;

end
