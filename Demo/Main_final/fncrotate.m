function [Xrtrain,Xrtest] = fncrotate(Xtrain,Xtest,W)
% Rota los datos y los Hankeliza
%

[~,covtr] = fncCov(Xtrain);

ntr = size(Xtrain,2);
nts = size(Xtest,2);
Xrtrain = cell(1,ntr);
Xrtest = cell(1,nts);

for tr = 1:ntr
%     [trace(covtr(:,:,tr))]
    Xtrain{tr} = Xtrain{tr}./sqrt(trace(covtr(:,:,tr)));
    Xrtrain{tr} = W*Xtrain{tr}';
end

clear covtr
[~,covtr] = fncCov(Xtest);
for tr = 1:nts
    Xtest{tr} = Xtest{tr}./sqrt(trace(covtr(:,:,tr)));
    Xrtest{tr} = W*Xtest{tr}';
end

