%  csp(Ca,Cb,cumvar,mode,K,Q)
 
 clear all; close all; clc
%%

% % Direccion de la base de datos
% SUBJECTS_DIR = 'D:\BCI';
SUBJECTS_DIR = 'D:\BCI';
% % Direccion del fold de las funciones
% addpath(genpath('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
% addpath(genpath('C:\Users\lfvelasquezm\Desktop\frank\functions'))
% addpath(genpath('C:\Users\frany\Dropbox\Event-related\Codes\TP\Matlab_wang\csp\CSP_fun\functions'));
addpath(genpath('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
%
% %% DataBase
% % BCIIII_4a_
% % BCICIV_2a_
%  GIGASCIENCE_
%
COHORT = 'BCICIV_2a_';
SUBJECTS = dir([SUBJECTS_DIR filesep '*' COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,:)';
%
%
% %% grilla de busqueda
param = linspace(0,0.9,100);

experiment_name = mfilename;

SS = 1:9;

%% paramaters definition
PPval = cell(numel(SS),1);
rho = [];
tstart = 0;
tend = 2;
% definir parametros de filter bank
f_low  = 4;
f_high = 40; 
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';
orden_filter = 5;
labels = [1 2];

%%
for s = SS
%     clearvars -except s SS rho experiment_name COHORT param SUBJECTS SUBJECTS_DIR Acc table PPval Rep tstart tend
    %
    load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'])
    y = y(:);
    ind = ismember(y,labels);
    y = y(ind);
    X = X(ind);
    X = cellfun(@(x) double(x) ,X,'UniformOutput',false);
    tic 
    ind = ismember(y,labels);
    %definitions
    acc=nan(5,numel(param));
    ks=nan(5,numel(param));
    Xcp = cell(5);
    sfeats = cell(5);
    Xa = cell(size(X,1),1);
    for fold = 1:5
        %         tic;
        tr_ind   = cv.training(fold); tr_ind = tr_ind(ind);
        ts_ind   = cv.test(fold); ts_ind = ts_ind(ind);
        
        Xa = cellfun(@(x) x(seg_start:seg_end,:),X,'UniformOutput',false);
        C = cell2mat(reshape(cellfun(@(x)(cov(x)/trace(cov(x))),Xa,'UniformOutput',false),[1 1 numel(Xa)]));
        W = csp_feats(C(:,:,tr_ind),y(tr_ind),'train','Q',floor(6/2));
        Xc = csp_feats(C,W,'test');
          
        %% LDA
        %
        C1 = cov(Xc(y==labels(1),:)); C2 = cov(Xc(y==labels(2),:));
        covTotal = C1 + C2;
        M1 = mean(Xc(y==labels(1),:)); M2 = mean(Xc(y==labels(2),:));
        SB = (M2-M1)'*(M2-M1);      
        [Ut,Dt] = eig(covTotal); %caution: the eigenvalues are initially in increasing order
        eigenvalues = diag(Dt);
        [eigenvalues,egIndex] = sort(eigenvalues, 'descend');
        Ut = Ut(:,egIndex);
        P = diag(sqrt(1./eigenvalues)) * Ut';
        %transforming covariance matrix of first class using P
        transformedCov1 =  P * SB * P';
        %EVD of the transformed covariance matrix
        [U1,D1] = eig(transformedCov1);
        eigenvalues = diag(D1);
        [eigenvalues,egIndex] = sort(eigenvalues, 'descend');
        U1 = U1(:, egIndex);
        W = U1' * P;
        W = W([1:Q end-Q+1:end],:);
        
       
        
            % Lasso
            target = mapminmax(y(tr_ind)')';
            B = lasso(Xc(tr_ind,:),target,'Lambda',param);
            selected_feats = abs(B)>eps;      
            sfeats{fold,u} = selected_feats;
            %
            for l=1:numel(param)
                Xcc = Xc(:,selected_feats(:,l));
                if size(Xcc,2)<2                    
                    continue
                end
                mdl = fitcdiscr(Xcc(tr_ind,:),y(tr_ind)); %LDA
                acc(fold,u,l) = mean(mdl.predict(Xcc(ts_ind,:))==reshape(y(ts_ind),[sum(ts_ind) 1]));
                %Confusion Matrix
                tar_pred = mdl.predict(Xcc(ts_ind,:)); %tar_pred(tar_pred==1)=0; tar_pred(tar_pred==2)=1;
                tar_true = reshape(y(ts_ind),[sum(ts_ind) 1]); %tar_true(tar_true==1)=0; tar_true(tar_true==2)=1;
                conM = confusionmat(tar_true,tar_pred);
                ks(fold,u,l) = kappa(conM);
                %plotconfusion(tar_true',tar_pred');
            end %lambda
            fprintf(['Threshold...' num2str(u) '...' num2str(toc) '\n'])
            %             [fold,u]
        %end % signficance
        fprintf(['Sujeto: ' num2str(s) ' Fold...' num2str(fold) '\n'])
%         toc
    end % folds 
    act = squeeze(mean(acc,1));
    [dato,indp] = max(act(:));
    actstd = squeeze(std(acc,1)); actstd = actstd(:); actstd = actstd(indp);
    [u_opt,l_opt]=ind2sub(size(act),indp);
    inds = SUBJECTS{s}; inds = inds(end-1:end); inds = str2double(inds);
    table(inds,:) = [threshold(u_opt),param(l_opt),dato*100,actstd*100];
%     table(inds,:) = param_subject(acc,threshold,param);
    toc
    %% Guardar resultados
%     save([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'results\' experiment_name 'acc.mat'],'acc','table');
   save(['D:\BCI' ...
        filesep SUBJECTS{s} filesep 'results\' experiment_name 'Results.mat'],'acc','table','Xcp','sfeats','ks');
    
   save(['C:\Users\lfvelasquezm\Dropbox\ERD\results_ERDfc_subjects\BCI' ...
        filesep SUBJECTS{s} filesep experiment_name 'Results.mat'],'acc','table','ks');
    fprintf([' ...acc: ' num2str(dato*100,'%02.1f') ' std: ' num2str(actstd*100,'%02.1f')...
        ' ...time: ' num2str(toc) '\n']);
end
