function [dd] = fncHank_(Xdr,e)
% dat es la matriz de la venta.
%dd = cell(1,numel(Xdr));
%for s = 1:numel(Xdr)
    sample = numel(Xdr);
%     dd{s} = cell(sample,1);
    for j = 1:sample
        x = Xdr{j};
        dd{j} = hankel(x(1:e),[x(e:end);x(1:e-1)]); % realiza con el vector la hankelizacion.
%         num = size(ddat,2);
        %dd{s}{j} = ddat(:,1:end-1);
    end
%end