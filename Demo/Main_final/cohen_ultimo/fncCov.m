function [covXr,covtr] = fncCov(Xdr,mode)
%input 
% Data reference
% mode = 1 - samples, 0- channels
%Output:

%%
% [nchan,nsamples] = size(Xdr{1});
% if nchan>nsamples; [nchan,~] = size(Xdr{1}'); end
ntr = length(Xdr);
%covtr = zeros(nchan,nchan,ntr);

for j = 1:ntr
    E = Xdr{j};
    if mode ~= 0; E = E'; end 
    covtr(:,:,j) = cov(E);
    covtr(:,:,j)= covtr(:,:,j)/trace(covtr(:,:,j));
end
covXr = mean(covtr(:,:,:),3);




