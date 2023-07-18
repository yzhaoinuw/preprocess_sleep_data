function [kappa] = cohenskappa(X,Y)

if nargin==2
    X = crosstab(X,Y);
end
    

% Recast X as proportions
X = X./sum(X(:));

% Observed Agreement
PRo = sum(diag(X)) / sum(X(:));

% Expected Agreement
PRe = sum(sum(X,1) .* sum(X,2)');

kappa = (PRo - PRe) / (1 - PRe);

if nargin < 1
    disp(['kappa = ',num2str(kappa,'% 4.3f')]);
end


