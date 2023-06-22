%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for different channels
% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SampEn] = Main_SE(series,ul1,up1,t3)
[N_trial,~] = size(series);
for trial = 1:N_trial
    %tic
    dat = squeeze(series(trial,:));
    a=1;
    for ti = ul1:up1
        %            x = SampEn2(2,0.45*std(dat(ti:ti+t3)),dat(ti:ti+t3));
        %            flag =  isinf(x);
        %            if flag == 1
        %                break
        %            end
        SampEn(trial,a) = SampEn2(2,0.45*std(dat(ti:ti+t3)),dat(ti:ti+t3));
        a = a +1;
    end
    %fprintf(['Channel: ' num2str(chan) num2str(toc) '\n'])
end