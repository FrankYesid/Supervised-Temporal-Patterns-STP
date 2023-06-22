sub=1:7;
ns=7;
np = 5; % Kfolds
[S,F] = ndgrid(sub,2:np);

avtrain = cell(size(S));
avtest = cell(size(S));
Yt = cell(1,ns);

for jj = 1:numel(S)
    s = S(jj);
    f = F(jj);
    Xd = X{s};
    % T-test
    fprintf('ttest \n')
    tt{s,f} = tes2(SampEnc1{s,f},SampEnc2{s,f});
    
    P_cmin = min(tt{s,f});  % Pc = min(f)-> Pcf
    P_cor = sort(P_cmin);   %ordena de menor a mayor
    
    % relaciona el numero del canal con el orden de los p-valores.
    for ch = 1:N_channel
        tmp = find(P_cmin == P_cor(ch));
        if length(tmp)>1
            tmp = tmp(1);           
            P_cmin(tmp) = 1000;
        end
        posch{s,f}(ch) = tmp;
    end
    time = toc 
    fprintf(['Sample Entropy done Sub: ' num2str(s) ' Fold:' num2str(f) '\n'])
    fprintf('Sub - %d. Fold - %d. Time: %f\n',s,f,time/60)
    save('ResultsSamEnBCIIV1_f2.mat','tt', 'posch' )
    %% CSP
    
   
    for ii = 6:59
        tic
        Xtrain_ = cellfun(@(x) x(:,posch{s,f}(1:ii)),Xtrain,'UniformOutput',false);
        Xtest_ = cellfun(@(x) x(:,posch{s,f}(1:ii)),Xtest,'UniformOutput',false);
        tini = 2.5; % (MI segment)
        trfin = 4.5;
        samples = (trfin-tini)*fs;
        Xtrain_ = fncCutdata_(Xtrain_,tini,trfin,fs);   % selecciona
        Xtest_ = fncCutdata_(Xtest_,tini,trfin,fs);   % selecciona
        
%         Xtrain_ = cellfun(@(Xa) Xa-mean(Xa),Xtrain_,'UniformOutput',false); % media cero    
%         Xtest_ = cellfun(@(Xa) Xa-mean(Xa),Xtest_,'UniformOutput',false); % media cero
         
        W = fncCSP(Xtrain_,ltr,3);
        
        %feature
        [Xftrain,Xftest] = fncfeatlog(Xtrain_,Xtest_,W');
        Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
        predtr = (predict(Mdl,Xftrain));
        predts = (predict(Mdl,Xftest));
        
        avtrain{jj}(ii) = sum(predtr == ltr)/numel(ltr);
        avtest{jj}(ii) = sum(predts == lts)/numel(lts);
        time = toc;
        fprintf('Suj: %d. Fold: %d. Relevence ch - %d.   Time: %f\n',s,f,ii,time)
        fprintf('Acc    : Training %f. Test %f\n', avtrain{jj}(ii)*100,avtest{jj}(ii)*100)
    end
end