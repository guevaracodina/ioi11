function out = ioi_cine_display_run(job)
out = [];
try
    %launch GUI
    ioi_cine_display_GUI;

catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end
end