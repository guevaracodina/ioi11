function IC = ioi_dfg_include_colors(OD,HbO,HbR,HbT,Flow)
include_flow      = cfg_menu;
include_flow.tag  = 'include_flow';
include_flow.name = 'Include flow';
include_flow.labels = {'Yes','No'};
include_flow.values = {1,0};
include_flow.val  = {Flow};
include_flow.help = {'Include flow.'}';

include_HbT      = cfg_menu;
include_HbT.tag  = 'include_HbT';
include_HbT.name = 'Include HbT';
include_HbT.labels = {'Yes','No'};
include_HbT.values = {1,0};
include_HbT.val  = {HbT};
include_HbT.help = {'Include HbT.'}';

include_OD      = cfg_menu;
include_OD.tag  = 'include_OD';
include_OD.name = 'Include optical intensity';
include_OD.labels = {'Yes','No'};
include_OD.values = {1,0};
include_OD.val  = {OD};
include_OD.help = {'If the optical intensity images (Green, Red, Yellow) have not been deleted'
    'previously, choose whether to generate movies for these colors.'}';

include_HbO      = cfg_menu;
include_HbO.tag  = 'include_HbO';
include_HbO.name = 'Include HbO';
include_HbO.labels = {'Yes','No'};
include_HbO.values = {1,0};
include_HbO.val  = {HbO};
include_HbO.help = {'Include HbO.'}';

include_HbR      = cfg_menu;
include_HbR.tag  = 'include_HbR';
include_HbR.name = 'Include HbR';
include_HbR.labels = {'Yes','No'};
include_HbR.values = {1,0};
include_HbR.val  = {HbR};
include_HbR.help = {'Include HbR.'}';

IC         = cfg_branch;
IC.name     = 'Include colors';
IC.tag    = 'IC';
IC.val     = {include_OD include_HbO include_HbR include_HbT include_flow}; 
IC.help    = {'Choose colors to include.'};

