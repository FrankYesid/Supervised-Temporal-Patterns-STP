function [covXr,covtr] = fncCov(Xdr)
%input 
% Data reference
%Output:

%%
nchan = size(Xdr{1},2);
ntr = length(Xdr);
covtr = zeros(nchan,nchan,ntr);
% covXr = zeros(nchan,nchan);
for j = 1:ntr
    E = Xdr{j};
    covtr(:,:,j) = cov(E);
    covtr(:,:,j)= covtr(:,:,j)/trace(covtr(:,:,j));
end
covXr = mean(covtr(:,:,:),3);




