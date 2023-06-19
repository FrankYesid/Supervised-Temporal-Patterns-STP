
%% MAIN IWCSP
% COHEN2017: Using spatiotemporal source separation to identify 
% prominent features in multichannel data without sinusoidal filters
%
% L.F. Velasquez-Martinez 
% -------------------------------------------------------------------------

clc, clear %all;
%
%
 
%% Dataset
% dt = 1,2,3:  1:BCICIV_2a
%
% for dt = 1:3
% dt = 1; % load data
% dataset = dt;

% dataset =0;
% load('G:\work_dic\Articulo2016\codes\BCICIIIIVa.mat');
% name = 'BCICIIIVa';
% fs = 100;
% sub = 1:5; %% Subject to analyze
% Xd_ = cell(1,length(sub));
% lb = cell(1,length(sub));
% labels = [1 2];
% for i = 1:length(sub)
%     ind{i} = ismember(y{i},labels);
%     Xd_{i} = X{i}(ind{i});
%     lb{i} = y{i}(ind{i});
% end
% % (MI segment)
% tini = 0.0001;
% tfin = 3.5;
% Xdr = fncCutdataf(Xd_,tini,tfin,fs,[8,30]);  % ventana de tiempo
% fprintf('Loaded data: %s\n',name);
% ns = size(Xd_,2);

%%
load('BCICIV_2a.mat');% raw data
name = 'BCICIV2a';
dataset = 1;
ns = size(X,1);
sub = 1:ns;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [1 2]; % clases
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i}); %belonging to labels
    lb{i} = y{i}(ind{i});
end
fprintf('Loaded data: %s\n',name);

%% Filter Band
clc;
a = [8 30];                                     % rango del filtro.
A = cell(9,1);                                  % celdas de la señal filtrada.
filt = 0;
% filtrado de cad sujeto.
if filt == 1
    for k = 1:9
        fprintf(['Sujeto: ' num2str(k) ' de 9' '\n'])
        Xfreq = fcnfiltband(X{k},fs,a,5);           % funcion Filtro.
        A{k} = Xfreq;
    end
    X = A;
end

%% Channels C3, C4
Xa_ = cell(1,ns);
XC3 = cell(2,ns);
Chan = 8; % channel
for i = 1:ns
   Xa_{i} = cellfun(@(X) X(:,Chan),Xd_{i},'UniformOutput',false);
end
Xd_ = Xa_; clear Xa_;


%% Selecting data for analysis
% (MI segment)
tini = 2;
tfin = 4;
samples = (tfin-tini)*fs;
Xdd = fncCutdata(Xd_,tini,tfin,fs);   % selecciona y filtra

% (Reference segment)
trini = 0;
trfin = 2;
samplesr = (trfin-trini)*fs;
Xdr = fncCutdata(Xd_,trini,trfin,fs);  % selecciona y filtra

%% hankelizar - time.
datd = fncHank(Xdd); % data
datr = fncHank(Xdr); %reference
Xdd = datd; clear datd;
Xdr = datr; clear datr;

%%
DatXc3 = fncHank2(Xd_); % data - time serial

%% Analysis
np = 5; % Kfolds
[S,F] = ndgrid(sub,1:np);
avtrain = zeros(size(S));
avtest = zeros(size(S));

% pruebas = numel(S);
%%
% 

for i = 1: numel(S)
    tic
    s = S(i);
    f = F(i); % fold
    
    Xd = Xdd{s};%Xdata
    Xr = Xdr{s};%Xreference
    
    tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
    ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
    
    Xtrd = Xd(tr_ind);
    Xtrr = Xr(tr_ind);
    
    Xtest = Xd(ts_ind);
    idtr = find(tr_ind==1);
    ltr = lb{s}(idtr)';
    idts = find(ts_ind==1);
    lts = lb{s}(idts)';
    
    %Etapa de diseño
    Wt = fncCSPcohen(Xtrd,Xtrr); %Finding temporal filters 
    %[Xd_,Xr_] = fncrotate(Xtrd,Xtrr,Ws); %Rotating data
    %[Xhd_,Xhr_] = fncHankel(Xd_,Xr_,500);
    %Wt = fncCSPcohen(Xhd_,Xhr_); %Finding temporal filters 
    
    %%
    Yt = cell(ns,1);
    for i = 1:ns
        load(['dd' num2str(i)  '.mat'])
        Yt{i} = cell(1,numel(dd{i}));
        for k = 1:numel(dd{i})
            Yt{i}{k} = cell(1,4);
            ja=1;
            for j = 1:500:1001
                Yt{i}{k}{ja} = Wt*dd{i}{k}(j:j+499,1:500);
                ja=ja+1;
            end
        end
        clear dd
    end
    
    %% Mean
    ym = cell(ns,1);
    for i = 1:ns
        Yt{i} = cell(1,numel(Yt{i}));
        for j = 1:numel(Yt{i})
            ym{i}{j} = mean(cell2mat(Yt{1}{1}'),1); % 1x1750
        end
%         ym{i} = mean(cell2mat(ym{i}));
    end
    
    %%
    %Aplicando rotaciones sobre Xtrain y Xtest
    [Xdtrain,Xdtest] = fncrotate(Xtrd,Xtest,Wt); %Rotating data
    [Xdftrain,Xdftest] = fncKernelfiltercohen(Xdtrain,Xdtest,Wt); %feature
    
    Mdl = fitcdiscr(Xdftrain,ltr); %%%%%%LDA%%%%%%
    predtr = (predict(Mdl,Xdftrain));
    predts = (predict(Mdl,Xdftest));
    
    if dataset == 1
        predtr = predtr';
        predts = predts';
    end
    avtrain(i) = sum(predtr == ltr)/numel(ltr);
    avtest(i) = sum(predts == lts)/numel(lts);
    tiempo = toc;
    fprintf('Exp - %d. Time: %f\n',i,tiempo)
    fprintf('Acc    : Training %f. Test %f\n', avtrain(i)*100,avtest(i)*100)
    %[i,k]
end
% % CSP
fprintf('Testing Results: %s \n', name)
for i = 1:numel(sub)
    acc_s(i,1) = mean(avtest(i,:,:))*100;
    acc_s(i,2) = std(avtest(i,:,:))*100;
    fprintf('Acc : %f. +- %f\n', acc_s(i,1),acc_s(i,2))
end

% 

%% 


%%
% %%%%%%%%%%%%%%%%%%%%%% Graphs

% save('results_CSP1.mat','av_s','std_s');
% figure;
% errorbar(1:9,av_s,std_s,'-s','LineStyle','none','MarkerSize',10,...
%     'MarkerEdgeColor','blue','MarkerFaceColor','blue'), ylim([0 1])
% figure
% bar(1:9,av_s)
% ylim([0 1])
%
% %% TRCSP (alpha)
% for i = 1:numel(sub)
%     for j = 1:numel(alph)
%     av_s(j,i) = mean(mean(squeeze(avtest(i,:,j))));
%     std_s(j,i) = std((squeeze(avtest(i,:,j))));
%     end
% end
% save('results_WTRCSP.mat','av_s','std_s','weightVector');
% figure;
% for j = 1:9
%     subplot(3,3,j)
%     fncGraph_acc(av_s(:,j)',std_s(:,j)',alph,0,1);
%     title(['S',num2str(j)]);
% end
% %% WTRCSP
% figure;
% for j = 1:9
%     plot(1:9,weightVector{j},'-s','MarkerSize',10,...
%      'MarkerEdgeColor','blue','MarkerFaceColor','blue'); ylim([0 1]);
%     title(['S',num2str(j)]);
% end
