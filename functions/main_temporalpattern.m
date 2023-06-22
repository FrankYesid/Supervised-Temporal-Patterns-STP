
%% MAIN IWCSP
% COHEN2017: Using spatiotemporal source separation to identify
% prominent features in multichannel data without sinusoidal filters
%
% L.F. Velasquez-Martinez
% F.Y. Zapatata-Castaño
% -------------------------------------------------------------------------

%% Databases
% load('D:\Luisa\Codes\BCICIV_2a.mat');               % raw data
% load('H:\ERD\Database\BCICIV_2a.mat');            % raw data
% name = 'BCICIV2a';
% load('D:\Archivos\Databasses_BCI\BCICIV_2a.mat')    % raw data
% load('H:\temporal_pattern\Load-data-struct\databasecastellano.mat');            % raw data
% name = 'Arg_cast';
% clc; clear;
load('BCICIV_1.mat')   % carga base de datos EEG. BCICIV_1
name  = 'BCIIV_1';
fprintf('Loaded data: %s\n',name);
ns = numel(X);
sub = 1:ns;
fs = 100;
Xddd = cell(1,ns);
lb = cell(1,ns);
N_channel = size(X{1}{1},2);
labels = [1 2];
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xddd{i} = X{i}(ind{i}); %belonging to labels
    lb{i} = y{i}(ind{i});
end

wwtt = cell(1,ns);
PP1 = cell(1,ns);


%% Analysis
np = 5; % Kfolds
[S,F] = ndgrid([1,2,6,7],1:np);

avtrain = cell(size(S));
avtest = cell(size(S));
jray = cell(size(S));
dkl = cell(size(S));

% SampEn1_ = cell(size(S));
% SampEn2_ = cell(size(S));
% tt = cell(size(S));
posch = cell(ns,1);

Yt = cell(1,ns);
%load('Results_TP_BCIIV1_f2.mat', 'Yt');
% load('Results_TP_BCIIV1.mat','Yt')

% load('ResultsSamEnBCIIV1_f2.mat','SampEn1_','SampEn2_','tt', 'posch' )
%_sub4_7
% load('tt.mat')
% P_cmintemp = cellfun(@(X) min(X),tt,'UniformOutput',false);  % Pc = min(f)-> Pcf
% P_cmin = cell(ns,1);
% for s = 1:ns
%     P_cmin{s} = cell2mat(P_cmintemp(s,:)');
% end
%
% P_cmin = cellfun(@(X) mean(X),P_cmin,'UniformOutput',false);
%
% P_cor = cellfun(@(x) sort(x),P_cmin,'UniformOutput',false);   %ordena de menor a mayor
%
% % relaciona el numero del canal con el orden de los p-valores.
% for s = 1:ns
%     for ch = 1:N_channel
%         tmp = find(P_cmin{s} == P_cor{s}(ch));
%         if length(tmp)>1
%             tmp = tmp(1);
%             P_cmin{s}(tmp) = 1000;
%         end
%         posch{s}(ch) = tmp;
%     end
% end
% 
% load('ResultSampEnBCIIV1_all_folds.mat')
% SampEn1_ = samEn1;
% % SampEn2_ = samEn2;
% for jj = 1:numel(S)
%     %     tic
%     s = S(jj);
%     f = F(jj);
%     Xd = X{s};
%     tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
%     ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
%     
%     Xtrain = Xd(tr_ind);
%     Xtest = Xd(ts_ind);
%     idtr = find(tr_ind==1);
%     ltr = lb{s}(idtr)';
%     idts = find(ts_ind==1);
%     lts = lb{s}(idts)';
%     
%     [~, N_channel] = size(Xtrain{1});
%         tic
%     
%         for labels = 1 :2 % clases
%             fprintf(['Sub: ' num2str(s) ' Fold' num2str(f) ' Clase: ' num2str(labels) '\n'])
%             Xd_ = Xtrain(ltr==labels); % belonging to labels
%             N_trials = numel(Xd_);
%             [sam, N_channel] = size(Xd_{1});
%     
%             Yttemp = zeros(N_trials,sam);
%             parfor c = 1:N_channel
%                 % solo canal seleccionado
%                 Xa_ = cellfun(@(N) N(:,c),Xd_,'UniformOutput',false);
%                 Xa_ = cellfun(@(Xa) Xa-mean(Xa),Xa_,'UniformOutput',false); % media cero
%     
%                 %% Selecting data for analysis
%                 tini = 2; % (MI segment)
%                 trfin = 4;
%                 samples = (trfin-tini)*fs;
%                 Xdc = fncCutdata_(Xa_,tini,trfin,fs);   % selecciona
%     
%                 trini = 0;  % (Reference segment)
%                 trfin = 2;
%                 samplesr = (trfin-trini)*fs;
%                 Xrc = fncCutdata_(Xa_,trini,trfin,fs);
%     
%                 w = 180; % hankelizar
%                 Xdd = fncHank_(Xdc,w); % data
%                 Xdr = fncHank_(Xrc,w); %reference
%     
%                 Wt = fncCSPcohen(Xdd,Xdr,2,0); % Finding temporal filters
%                 Wtemp(c,:) = Wt;
%                 Ptemp(c,:) = fncFFT(Wt,fs);
%                 % Yt = Wt*Xc
%                 Yt_conv = cellfun(@(X) conv(X,Wt','same'),Xa_,'UniformOutput',false);
%                 Yttemp(:,:,c) = cell2mat(cellfun(@(X) envelope(X)',Yt_conv,'UniformOutput',false)');
%     
%                 %             fprintf(['Sub: ' num2str(s) ' Fold: ' num2str(f) ' Clase: ' num2str(labels) ' Canal: ' num2str(c) '\n'])
%             end
%             wwtt{s,f}{labels} = Wtemp;
%             PP1{s,f}{labels} = Ptemp;
%             Yt{s,f}{labels} = Yttemp;
%         end
%         time = toc;
%         fprintf(['Temporal patterns done Sub: ' num2str(s) ' Fold: ' num2str(f) '\n'])
%         fprintf('Exp - %d. Time: %f\n',jj,time)
%         save('Results_TP_BCIIV1_all.mat', 'Yt','PP1','wwtt');
%     
    %% Relevance
    load('Results_TP_BCIIV1_all.mat', 'wwtt')
    for s = [1,2,6,7]
        for fold = 1:5
            for cla = 1:2
                tic 
                for ch = 1:59
                    temp = wwtt{s,fold}{cla}(ch,:);
                    SampEn_{s,cla}(fold,ch) = SampEn2(2,0.45*std(temp),temp);
                end
                toc 
            end
        end
    end
    %%
    for s = [1,2,6,7]
        figure
        for fold = 1:5
            subplot(1,5,fold)
            for cla = 1:2
                boxplot([SampEn_{s,cla}(fold,:)',SampEn_{s,cla}(fold,:)'],'Notch','on','Labels',{'class 1','class 2'})
            end
            title(['Subject ' num2str(s) ' - Fold ' num2str(fold)])
        end
        ylim([0 1])
    end
    %%
    % Samp Entropy
    % Calculate SampEn using a sliding time window of width 2s from 2.5 to 3.5s for each class.
    % Analysis window
%         ul1 = 250;                                        % ul inicial.
%         up1 = 350;                                        % up final.
%         w = 200; % seg
    %     tic
    %     fprintf('SampEn1')
    %     SampEn1_{s,f} = Main_SE(Yt{s,f}{1},ul1,up1,w);
    %     fprintf('SampEn2')
    %     SampEn2_{s,f} = Main_SE(Yt{s,f}{2},ul1,up1,w);              % ERD class 2
    %
    %     % T-test
    %     fprintf('ttest')
%%
  
%     tt{s,f} = tes2(SampEn1_{s,f},SampEn2_{s,f});
    
%     P_cmin = min(tt{s,f});  % Pc = min(f)-> Pcf
%     P_cor = sort(P_cmin);   %ordena de menor a mayor
%     
%     % relaciona el numero del canal con el orden de los p-valores.
%     for ch = 1:N_channel
%         tmp = find(P_cmin == P_cor(ch));
%         if length(tmp)>1
%             tmp = tmp(1);
%             P_cmin(tmp) = 1000;
%         end
%         posch{s,f}(ch) = tmp;
%     end
%     time = toc;
%     fprintf(['Sample Entropy done Sub: ' num2str(s) ' Fold:' num2str(f) '\n'])
%     fprintf('Sub - %d. Fold - %d. Time: %f\n',s,f,time/60)
%     %     save('ResultsSamEnBCIIV1_22.mat','SampEn1_','SampEn2_','tt', 'posch' )
    
    
    
    %% CSP
    
    %
    %     for ii = 6:59
    %         tic
    %         Xtrain_ = cellfun(@(x) x(:,posch{s}(1:ii)),Xtrain,'UniformOutput',false);
    %         Xtest_ = cellfun(@(x) x(:,posch{s}(1:ii)),Xtest,'UniformOutput',false);
    %         tini = 2.5; % (MI segment)
    %         trfin = 4.5;
    %         samples = (trfin-tini)*fs;
    %         Xtrain_ = fncCutdata_(Xtrain_,tini,trfin,fs);   % selecciona
    %         Xtest_ = fncCutdata_(Xtest_,tini,trfin,fs);   % selecciona
    %
    %                 Xtrain_ = cellfun(@(Xa) Xa-mean(Xa),Xtrain_,'UniformOutput',false); % media cero
    %                 Xtest_ = cellfun(@(Xa) Xa-mean(Xa),Xtest_,'UniformOutput',false); % media cero
    %
    %         [W,Cov] = fncCSP(Xtrain_,ltr,3);
    %
    %                 parfor ii = 1:numel(F) %frequencies and time
    %         Xd_ = X{ii};
    %         [W,Cov] = fncCSP(Xd_,ytr,3); %------------ W rotation CSP-based
    %         ------------ Rayleight quotient using the first component
    %         jtemp(ii) = trace(W*Cov{1}*W')/trace(W*Cov{2}*W');
    %
    %         feature
    %         [Xftrain,Xftest,covXtest] = fncfeatlog(Xtrain_,Xtest_,W');
    %         Xf1 = Xftest(lts==1,:);
    %         Xf2 = Xftest(lts==2,:);
    %         Cov{1} = covXtest(:,:,lts==1); Cov{1}=mean(Cov{1},3);
    %         Cov{2} = covXtest(:,:,lts==2);  Cov{2}=mean(Cov{2},3);
    %
    %         Div
    %         k = 6;
    %         for ff = 1 :6
    %         tmp1 = trace((cov(Xf1))\cov((Xf2)));
    %         tmp2 = (mean(Xf1)-mean(Xf2))*inv(cov(Xf2))*(mean(Xf1)-mean(Xf2))';
    %         tmp3 = log(det(cov(Xf1))/det(cov(Xf2)));
    %         tdkl = 0.5*(tmp1+tmp2-k+tmp3);
    %         end
    %         dkl{jj}(ii) = sum(tdkl);
    %         dkl{jj}(ii) = tdkl;
    %
    %         Rayleight
    %         j(ii) = (W(1,:)*Cov{1}*W(1,:)')/(W(1,:)*Cov{2}*W(1,:)');
    %         jtemp(ii) = (W(:,1)'*Cov{1}*W(:,1))/(W(:,1)'*Cov{2}*W(:,1));
    %         trace(W*Cov{1}*W')/trace(W*Cov{2}*W');
    %         Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
    %         predtr = (predict(Mdl,Xftrain));
    %         predts = (predict(Mdl,Xftest));
    %
    %         avtrain{jj}(ii) = sum(predtr == ltr)/numel(ltr);
    %         avtest{jj}(ii) = sum(predts == lts)/numel(lts);
    %         time = toc;
    %         fprintf('Suj: %d. Fold: %d. Relevence ch - %d.   Time: %f\n',s,f,ii,time)
    %         fprintf('Acc    : Training %f. Test %f\n', avtrain{jj}(ii)*100,avtest{jj}(ii)*100)
    %     end
    %     jray{jj} = jtemp;
    %     save('ResultsAccBCIIV1_f2.mat','avtest','dkl','jray')
% end

% load('ResultsAccBCIIV1_f2.mat')
set(0,'DefaultFigureWindowStyle','docked')
for s = [1,2,6,7] %[3,4]
    for i =1:5; Dkl{s}(i,:) = dkl{s,i}/max(dkl{s,i}); end
    for i =1:5; j_{s}(i,:) = jray{s,i}/max(jray{s,i}); end
    for i =1:5; ac{s}(i,:) = avtest{s,i}; end
    figure
    % a = a+1;
    subplot(3,1,1)
    errorbar(mean(ac{s}),std(ac{s})); title(['Acc' ' subject: ' num2str(s)])
    subplot(3,1,2)
    errorbar(mean(Dkl{s}),std(Dkl{s})); title(['Div' ' subject: ' num2str(s)])
    ylim([0 1])
    subplot(3,1,3)
    errorbar(mean(j_{s}),std(j_{s})); title(['Rayleight' ' subject: ' num2str(s)])
end

% for s = [3,4]
% for i =1:4; j_{s}(i,:) = jray{s,i}; end
% figure
% errorbar(mean(j_{s}),std(j_{s}))
% end


xpru = (trace((W*Cov{2}*W')'*(W*Cov{1}*W')'))/(sqrt(trace(W*Cov{2}*W'*W*Cov{2}*W'))*sqrt(trace(W*Cov{1}*W'*W*Cov{1}*W')));

norm(diag(W*Cov{1}*W')-diag(W*Cov{1}*W'))

set(0,'DefaultFigureWindowStyle','docked')
for s = [1,2,6,7]
    figure
    % cl = gray(59);
    % cl = flip(cl);
    for i =1:5; ac{s}(i,:) = avtest{s,i}; end
    [~,chan] = max(mean(ac{s}));
    posc = posch{s}(1:chan);
    % for ic =1:59
    scatter3(pos(posc,1),pos(posc,2),pos(posc,3),300,'filled')%,'MarkerFaceColor',[cl(ic,:)]);
    %     hold on
    % end
    view(0,90)
    title(['Subject ' num2str(s)])
    text(pos(:,1),pos(:,2),pos(:,3),electrodes);
end




