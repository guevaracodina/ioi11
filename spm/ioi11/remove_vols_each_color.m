function remove_vols_each_color(IOI,str_color,f1,s1)
c1 = find(IOI.color.eng==str_color);
if length(IOI.sess_res{s1}.fname)>=c1
    fname_list = IOI.sess_res{s1}.fname{c1};
    if ~isempty(fname_list)
        fname = fname_list{f1};
        try
            delete(fname);
        catch
            disp(['Could not delete ' fname]);
        end
    end
end

