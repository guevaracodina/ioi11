function SCKS=ioi_SCKS_set_priors(SCKS)
%load priors
switch SCKS.PS.PhysioModel_Choice
    case 0 %Buxton-Friston
        SCKS = ioi_SCKS_priors_Buxton(SCKS,5);
%     case 1 %Zheng-Mayhew
%         SCKS = ioi_HDM_priors_ZM(SCKS,9); 
%     case 2 %Huppert1
%         SCKS = ioi_HDM_priors_Hu(SCKS,9);
%     otherwise
end