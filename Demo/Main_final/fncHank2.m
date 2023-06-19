function [dd] = fncHank2(data)
% Hankel.
for s = 1:numel(data)
    dd = cell(1,numel(data));
    tic
    sample = numel(data{s});
    dd{s} = cell(sample,1);
    for j = 1:sample
        drt = hankel(data{s}{j}',[data{s}{j}(end)' data{s}{j}']);
        dd{s}{j} = drt(1:1750:1:1750);
    end
    time = toc;
    fprintf(['Subject: ' num2str(s) ' de 9 - ' 'time: ' num2str(time) '\n'])
    save(['dd' num2str(s) '.mat'],'dd','-v7.3')
    clear dd;
end
dd = 0;