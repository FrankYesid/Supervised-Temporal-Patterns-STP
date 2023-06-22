function [Xftrain,Xftest,covXtest] = fncfeatlog(Xtrain,Xtest,W)
% [lotte2011]
Xftrain = zeros(size(Xtrain,2),size(W,2));
Xftest = zeros(size(Xtest,2),size(W,2));

for tr = 1:size(Xtrain,2)
%     proyTr = W'*Xtrain{tr}';
%     %generating the features as the log variance of the projected signals
%     variances = var(proyTr,0,2);
%     for f=1:length(variances)
%         Xftrain(tr,f) = log(variances(f));
%     end
    C = cov(Xtrain{tr});
    C = C/trace(C);
    tmp = W'*C*W;
%     Xftrain(tr,:) = log(tmp);
    Xftrain(tr,:) = log(diag(tmp)/trace(tmp));
    covXtrain(:,:,tr) = tmp;
end

for tr = 1:size(Xtest,2)
%     proyTr = W'*Xtest{tr}';
%     %generating the features as the log variance of the projected signals
%     variances = var(proyTr,0,2);
%     for f=1:length(variances)
%         Xftest(tr,f) = log(variances(f));
%     end
    C = cov(Xtest{tr});
    C = C/trace(C);
    tmp = W'*C*W;
    Xftest(tr,:) = log(diag(tmp)/trace(tmp));
    covXtest(:,:,tr) = tmp;
%     Xftest(tr,:) = log(tmp);
end