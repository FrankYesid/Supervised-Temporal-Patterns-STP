function [W] = fncCSPcohen(Xd,Xdr)
%Input:
%Xd : X data
%Xr : X reference

%Output:
%W: the learnt CSP filters (a [Nc*Nc] matrix with the filters as rows)

[covXd,~]= fncCov(Xd);               % Data covariance
[covXr,~] = fncCov(Xdr);             % Reference covariance

% covXd = covXd./max(max(covXd));
% covXr = covXr./max(max(covXr));

% covXd = (log(abs(covXd)));
% covXr = (log(abs(covXr)));

%whitening transform of total covariance matrix
[Ut,Dt] = eig(covXd); %caution: the eigenvalues are initially in increasing order
eigenvalues = diag(Dt);
[eigenvalues,egIndex] = sort(eigenvalues, 'descend');
Ut = Ut(:,egIndex);
P = diag(sqrt(1./eigenvalues)) * Ut';

%transforming covariance matrix of first class using P
transformedCov1 =  P * covXr * P';

%EVD of the transformed covariance matrix
[U1,D1] = eig(transformedCov1);
eigenvalues = diag(D1);
[eigenvalues,egIndex] = sort(eigenvalues, 'descend');
% U = U1;

U1 = U1(:, egIndex([1]));

W = U1' * P;