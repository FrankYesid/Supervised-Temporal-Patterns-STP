clc; clear all;
load('BCICIV_1.mat')   % carga base de datos EEG. BCICIV_1
name  = 'BCIIV_1';
fprintf('Loaded data: %s\n',name);
ns = numel(X);
sub = 1:ns;
%ns = length(sub);
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


load('tt_final.mat')
P_cmintemp = cellfun(@(X) min(X),tt,'UniformOutput',false);  % Pc = min(f)-> Pcf
P_cmin = cell(ns,1);
for s = [1,2,6,7]%1:ns
    P_cmin{s} = cell2mat(P_cmintemp(s,:)');
end

P_cmin = cellfun(@(X) mean(X),P_cmin,'UniformOutput',false);

P_cor = cellfun(@(x) sort(x),P_cmin,'UniformOutput',false);   %ordena de menor a mayor

% relaciona el numero del canal con el orden de los p-valores.
for s = [1,2,6,7]
    for ch = 1:N_channel
        tmp = find(P_cmin{s} == P_cor{s}(ch));
        if length(tmp)>1
            tmp = tmp(1);
            P_cmin{s}(tmp) = 1000;
        end
        posch{s}(ch) = tmp;
    end
end
%% ray

%% Analysis
np = 5; % Kfolds
[S,F] = ndgrid([2],1:np);
%[1,2,6,7]
avtrain = cell(size(S));
avtest = cell(size(S));
% jray = cell(size(S));
CR = cell(size(S));
tini = 0.02; % 2.5;
tfin = 8; %4.5;
window = floor(fs*1);%0.125
overlap = 0.90;
lags = window+1:window-round(overlap*window):tfin*fs-tini*fs;
avtraintemp = zeros(1,length(lags));
avtesttemp = zeros(1,length(lags));
vchans = [6,18,59;6,35,59;6,19,59;6,19,59]; 


for jj = 1:numel(S)
    %     tic
    s = S(jj);
    f = F(jj);
    Xd = X{s};
    tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
    ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
    
    Xtrain = Xd(tr_ind);
    Xtest = Xd(ts_ind);
    idtr = find(tr_ind==1);
    ltr = lb{s}(idtr)';
    idts = find(ts_ind==1);
    lts = lb{s}(idts)';
    [~, N_channel] = size(Xtrain{1});
    
    %% CSP
    if s == 1; vchans = [6,18,59];        
    elseif s==2; vchans = [6,35,59];
    elseif s==6; vchans = [6,19,59];
    elseif s==7; vchans = [6,19,59]; 
    end
    for ii = vchans
        tic
        Xtrain_ = cellfun(@(x) x(:,posch{s}(1:ii)),Xtrain,'UniformOutput',false);
        Xtest_ = cellfun(@(x) x(:,posch{s}(1:ii)),Xtest,'UniformOutput',false);
        
        for t_ = 1:length(lags)
            seg = lags(t_);
            Xtrain__ = fncCutdataf(Xtrain_,(seg-window)/fs,(seg-1)/fs,fs,[8,30]);   % selecciona
            Xtest__ = fncCutdataf(Xtest_,(seg-window)/fs,(seg-1)/fs,fs,[8,30]);   % selecciona
            %         Xtrain_ = cellfun(@(Xa) Xa-mean(Xa),Xtrain_,'UniformOutput',false); % media cero
            %         Xtest_ = cellfun(@(Xa) Xa-mean(Xa),Xtest_,'UniformOutput',false); % media cero
            
            [W,Cov] = fncCSP(Xtrain__,ltr,3);
            
            %         parfor ii = 1:numel(F) %frequencies and time
            %         Xd_ = X{ii};
            %         [W,Cov] = fncCSP(Xd_,ytr,3); %------------ W rotation CSP-based
            %------------ Rayleight quotient using the first component
            %         jtemp(ii) = trace(W*Cov{1}*W')/trace(W*Cov{2}*W');
            
            %feature
            [Xftrain,Xftest,covXtest] = fncfeatlog(Xtrain__,Xtest__,W');
            Xf1 = Xftest(lts==1,:);
            Xf2 = Xftest(lts==2,:);
            Cov{1} = covXtest(:,:,lts==1); Cov{1}=mean(Cov{1},3);
            Cov{2} = covXtest(:,:,lts==2);  Cov{2}=mean(Cov{2},3);
            
            %Rayleight
            %         j(ii) = (W(1,:)*Cov{1}*W(1,:)')/(W(1,:)*Cov{2}*W(1,:)');
            jtemp(t_) = (W(:,1)'*Cov{1}*W(:,1))/(W(:,1)'*Cov{2}*W(:,1));
            Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
            predtr = (predict(Mdl,Xftrain));
            predts = (predict(Mdl,Xftest));
            
            avtraintemp(t_) = sum(predtr == ltr)/numel(ltr);
            avtesttemp(t_) = sum(predts == lts)/numel(lts);
            time = toc;
            fprintf('Suj: %d. Fold: %d. Relevence ch - %d.   Time: %f\n',s,f,ii,time)
            fprintf('Acc    : Training %f. Test %f.  Ventana. %d\n', avtraintemp(t_)*100,avtraintemp(t_)*100,t_)
        end
        
        avtrain{jj}{ii} = avtraintemp;
        avtest{jj}{ii} = avtraintemp;
        jray{ii} = jtemp;    
    end
    CR{jj} = jray; % fold sujetos - ventanas
    save('ResultsCRBCIIV1_s2.mat','CR')
end

clear jRcm
for s = 1
%     if s == 1; vchans = [6,18,59];        
    %if s==1; vchans = [6,35,59];
    if s==6; vchans = [6,19,59];
    elseif s==7; vchans = [6,19,59]; 
    end
    cc=1;
    for a = vchans
        for fl = 1:5
        jfc{s,cc}(fl,:) = CR{s,fl}{a}./max(CR{s,fl}{a});
        end
        jRcm{s,cc} = mean(jfc{s,cc});
        cc = cc+1;
    end
end

%%
set(0,'DefaultFigureWindowStyle','docked')
for s = 1
    figure
    if s == 1; vchans = [6,18,59]; s_=1;
    elseif s==2; vchans = [6,35,59]; s_=2;
    elseif s==3; vchans = [6,19,59]; s_=6;
    elseif s==4; vchans = [6,19,59]; s_=7;
    end
    for a = 1:3
        plot(lags/fs,jRcm{s,a}); hold on
    end
    
%     set(gca,'YScale','log')
set(gca,'YScale','log')
ylim([0 1])
    legend({'6 Channels', [num2str(vchans(2)) ' Channels'],'59 Channels'},'FontSize',9,'Orientation','horizontal','Interpreter','latex')
    legend('boxoff')
    title(['Subject' num2str(s_)])
end

