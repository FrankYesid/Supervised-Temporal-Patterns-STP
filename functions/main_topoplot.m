% load('mycamp')

% load position.mat
load HeadModel.mat
load('electrodesBCICIV1.mat')
load('final.mat')
rel = 1:22; % para colocar mas importancia al electrodo.
sub = [18,34,19,19];
a=1;
for s = [1,2,6,7]
%     subplot(1,4,a)
    figure(s)
    set(0,'DefaultFigureWindowStyle','docked')
    sel = posch{s}(1:sub(a)); % que canales selecciono.
    M1.xy = pos(:,1:2);
    M1.lab = electrodes;
    MyTopo_fun(rel,sel,M1.xy,M1.lab,[min(rel) max(rel)],0,0)
    caxis([-1,1])
    % colormap(gca,cmap);
    title(['Sujeto ' num2str(s)])
    axis square
    axis off
    a=a+1;
end