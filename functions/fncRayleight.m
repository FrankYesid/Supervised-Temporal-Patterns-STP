function [j] = fncRayleight(XdT,tr_ind,fs,twin,tseg,frl,bw,ytr,vt,vf)
% function: Computing the Rayleight measure discrimation indicator 
%
%
% L.F. Velasquez-Martinez
%------------ Computing Rayleight Quotient
[F,T] = ndgrid(frl,tseg);
X = cell(size(F));

Xd = XdT(tr_ind); % selecting just training trials

parfor ii = 1:numel(F) %------------ TIME-FRECUENCY GRID
    %         tic
    ta = T(ii);  %------------ time segment selection
    ti = ta/fs;
    tf = (ta+twin)/fs;
    
    fr = F(ii);  %------------ frequency band selection
    vfr = [fr,fr+bw];
    
    Xc = fncCutdataSub(Xdf,ti,tf,fs);
    
    X{ii} = Xc;
end
clear F T i

[F,~] = ndgrid(vf,vt); % freq and time
j = zeros(size(F));

parfor ii = 1:numel(F) %frequencies and time
    Xd_ = X{ii};
    [W,Cov] = fncCSP(Xd_,ytr,3); %------------ W rotation CSP-based
    %------------ Rayleight quotient using the first component
    j(ii) = (W(1,:)*Cov{1}*W(1,:)')/(W(1,:)*Cov{2}*W(1,:)');
end
