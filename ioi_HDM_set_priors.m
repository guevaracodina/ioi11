function [M]=ioi_HDM_set_priors(M)
%load priors
switch M.Model_Choice
    case 0 %Buxton-Friston
        M = ioi_HDM_priors_Buxton(M,5);
    case 1 %Zheng-Mayhew
        M = ioi_HDM_priors_ZM(M,9);   %MODIFIER
    case {2,7}%Huppert1
        M = ioi_HDM_priors_Hu(M,9);
%     case 3 %decay Z-M
%         M= ioi_HDM_priors_ZM(M,5);   %MO
%     case {4,5,6,8} %decay Hu
%         M= ioi_HDM_priors_Buxton(M,5);   %MO  
    otherwise
        %error('12')
end
%M.mAff=sum(strncmp(M.name,'eff',3)); % number of efficiency parameters to plot

% %reload pE and pC from calculated covariance
% if  M.priors(1)==2
%     M =doCovJob([],M,[]);
% elseif   M.priors(1)==3 % priors individual
%      M =doCovJob([],M,[]); % take covariance from LOU
%      load([M.rep 'donneesextraites.mat' ],'SPM2');    
%      M.var=SPM2(M.priors(2)).data.M.varEst;
% elseif   M.priors(1)==5 % priors individual
%      M =doCovJob([],M,[]); % take covariance from LOU     
%      load([strrep(M.rep,'lyse_1','lyse_3') 'donneesextraites.mat' ],'SPM2');    
%      M.var=SPM2(M.priors(2)).data.M.varEst;
% end



    
