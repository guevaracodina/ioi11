function flow_lambda_norm = ioi_flow_norm_cfg(choose_on)
% Graphical interface configuration options to normalize decorrelation velocity.
% 
% Reference:
% J. D. Briers and S. Webster, “Quasi real-time digital version of
% single-exposure speckle photography for full-field monitoring of velocity or
% flow fields,” Optics Communications, vol. 116, no. 1–3, pp. 36–42, Apr. 1995.
% 
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

IR_laser_lambda             = cfg_entry; 
IR_laser_lambda.name        = 'Laser Wavelength';
IR_laser_lambda.tag         = 'IR_laser_lambda';       
IR_laser_lambda.strtype     = 'r';
IR_laser_lambda.num         = [1 1];     
IR_laser_lambda.val         = {};
IR_laser_lambda.def         = @(val)ioi_get_defaults('flow1.lambda', val{:});
IR_laser_lambda.help        = {'Enter wavelength of laser in m (default value = 785e-9 m).'};

flow_lambda_norm_On         = cfg_branch;
flow_lambda_norm_On.tag     = 'flow_lambda_norm_On';
flow_lambda_norm_On.name    = 'Flow normalization ON';
flow_lambda_norm_On.val     = {IR_laser_lambda}; 
flow_lambda_norm_On.help    = {'Decorrelation velocity normalized by wavelength [Vc=lambda/(2*pi*Tc)]'};

flow_lambda_norm_Off        = cfg_branch;
flow_lambda_norm_Off.tag    = 'flow_lambda_norm_Off';
flow_lambda_norm_Off.name   = 'Flow normalization OFF';
flow_lambda_norm_Off.val    = {}; 
flow_lambda_norm_Off.help   = {'Decorrelation velocity not normalized [Vc=1/Tc]'};

flow_lambda_norm            = cfg_choice;
flow_lambda_norm.tag        = 'flow_lambda_norm';
flow_lambda_norm.name       = 'Flow Normalization';

flow_lambda_norm.values     = {flow_lambda_norm_On flow_lambda_norm_Off};
if choose_on
    flow_lambda_norm.val    = {flow_lambda_norm_On};
else
    flow_lambda_norm.val    = {flow_lambda_norm_Off};
end
flow_lambda_norm.help       = {'Choose whether to normalize the decorrelation velocity by the wavelength of the laser in speckle imaging'}';

% EOF
