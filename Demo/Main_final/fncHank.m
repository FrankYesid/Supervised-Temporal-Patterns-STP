function [dd] = fncHank(Xdr)
% dat es la matriz de la venta.
dd = cell(1,numel(Xdr));
for s = 1:numel(Xdr)
    sample = numel(Xdr{s});
    dd{s} = cell(sample,1);
    for j = 1:sample
        ddat = hankel(Xdr{s}{j}',[Xdr{s}{j}(end)' Xdr{s}{j}']);
        dd{s}{j} = ddat(1:500,1:500);
    end
end