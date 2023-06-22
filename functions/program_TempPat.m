clear; clc
load('final.mat','posch')                   % canales seleccionados.
load('Results_TP_BCIIV1_all.mat','wwtt')    % patrones temporales.
load('electrodesBCICIV1.mat', 'electrodes') % Nombre de los electrodos.
load('ResultSampEnBCIIV1_all_folds.mat')    % SampEn Clase 1 y 2.

fs = 100; % frecuencia de muestreo.

% Patron temporal
for s = [1,2,6,7]
    if     s == 1 Nchannels = [6,18,59];
    elseif s == 2 Nchannels = [6,34,59];
    elseif s == 6 Nchannels = [6,19,59];
    elseif s == 7 Nchannels = [6,19,59];
    end
    a = 1;
    for n_ch = Nchannels
        temp = wwtt(s,1:5);
        mm = cellfun(@(X) cellfun(@(x) x(posch{s}(1:n_ch),:),X,'UniformOutput',false),temp,'UniformOutput',false);
        mat = cellfun(@(X) cell2mat(X'),mm,'UniformOutput',false);
        mat = cell2mat(mat');
        tp{s,a} = mean(mat);
        a = a+1;
    end
end

%% plot TP
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
a=1;
for s = [1,2,6,7]                                                            % graficador en 3D todos los canales
    figure
%     subplot(2,4,a)
    hold on
    for ven = 1:numel(Nchannels)
        hold on
        plot(1/100:1/100:180/100,tp{s,ven},'linewidth',2)                    %% sin normalizar
    end
    yticks([])
    xticks([0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8])
    xticklabels([0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8])
    if     s == 1 legend({'6 Channels','18 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    elseif s == 2 legend({'6 Channels','34 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    elseif s == 6 legend({'6 Channels','19 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    elseif s == 7 legend({'6 Channels','19 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    end
    legend('boxoff')
    xlim([0 1.8])
    xlabel('Tiempo (s)')
    title(['Subject ' num2str(s)],'Interpreter','latex')
    drawnow
    hold off
    a = a+1;
end

%% Frecuencia
for s = [1,2,6,7]
    tem = cellfun(@(X) fncFFT(X,fs),tp(s,:),'UniformOutput',false);
    freq(s,:) = tem;
end

%% plot frecuencias - wt
set(0,'DefaultFigureWindowStyle','docked') % grafica cada figura en una misma pantalla.
a = 5;
for s = [1,2,6,7]   % graficador en 3D todos los canales
    figure
    L = 180; Fs = 100;
    f = Fs*(0:(L/2))/L;
%     subplot(2,4,a)
    hold on
    for ven = 1:numel(Nchannels)
        plot(f,freq{s,ven})%,'linewidth',2);%
    end
%     set(gca,'XScale','log')
    
%     yticks([])
    if     s == 1 legend({'6 Channels','18 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    elseif s == 2 legend({'6 Channels','34 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    elseif s == 6 legend({'6 Channels','19 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    elseif s == 7 legend({'6 Channels','19 Channels','59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    end
    legend('boxoff')
    xlim([8 30])
    ylim([10^(-4.5) 10^(0)])
    set(gca,'YScale','log')
    xlabel('Frequency (Hz)','Interpreter','latex')
    title(['Subject ' num2str(s)],'Interpreter','latex')
    drawnow
    hold off
    a = a+1;
end

%% 


