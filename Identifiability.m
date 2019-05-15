%% CHEME 5999 Pset 3 -> Identifiability

% this code determines which parameters are identifiable based on SSM 

function [List,pset] = Identifiability(SSM)

P_i_name = char('Km1','Km2','kcat','n'); % The list of all parameters in order

[r,c] = size(SSM); % r=time steps, c=4 
X =[];
pset =[];

% McAuley procedure doi:10.1081/PRE-120024426
R = SSM; % the first time, use SSM to do column sumComputeIdentifiability

for j = 1:c
     for i = 1:c
         M(i)= R(:,i)'*R(:,i); % the square sum of each column
     end
     
     [a,pos]= max(M)          % finds the indices of the maximum values of M, and returns them in output vector pos. 
     
     if a>1.0e-08                 % this is the tolerance
         X = [X SSM(:,pos)];      % the column that has the largest SS magnitude
         pset = [pset; pos];      % give the index of the parameter
         Shat = X*inv(X'*X)*X'*SSM; % Find the prediction SSM
         R = SSM-Shat;            % residual mtrx, now the residual matrix is the new mtrx that we find the next identifiable parameter, return to the top of the J loop  
     end
end

pset = unique(pset);                % identifiable parameters (w/o repeated values!) 

c = size(pset); % number of identifiable parameters 

for i = 1:c
    k = pset(i);
    List(i,:) = P_i_name(k,:);
end
 


end  