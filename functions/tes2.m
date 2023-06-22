function [tt] = tes2(A,B)
[trials,tau] = size(A);
for tr = 1:trials
    [h p] = ttest2(A(tr,:),B(tr,:));
    tl(1,1) = h;
    tl(1,2) = p;
    tt(tr,:) = tl;
end