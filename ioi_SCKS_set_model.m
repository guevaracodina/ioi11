function SCKS = ioi_SCKS_set_model(SCKS)
% Physiological Model choice
switch SCKS.PS.PhysioModel_Choice
    case 0 %Buxton-Friston
        SCKS.f     = 'ioi_fx_SCKS';
        SCKS.g     = 'ioi_gx_SCKS';
        SCKS.x     = zeros(4,1);       
    case 1 %Zheng-Mayhew
        SCKS.f     = 'ioi_HDM_fx_ZM';  
        SCKS.g     = 'ioi_HDM_gx';  
        SCKS.x     = zeros(5,1);    
    case 2 %Huppert1
        SCKS.f     = 'ioi_HDM_fx_Hu2';  
        SCKS.g     = 'ioi_HDM_gx';  
        SCKS.x     = zeros(7,1);           
end
SCKS.n     = length(SCKS.x);
%Number of inputs of direct model
SCKS.m=1;

