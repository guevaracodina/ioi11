function ioi_write_excel_core(writeET,v,sheet,tit,dir_fig)
if writeET == 1 || writeET == 3
    File = fullfile(dir_fig,[tit '.xls']); %.xls or .xlsx
    try
        Excel = actxserver('Excel.Application');
        if ~exist(File,'file')
            ExcelWorkbook = Excel.workbooks.Add;
            ExcelWorkbook.SaveAs(File,1);
            ExcelWorkbook.Close(false);
        end
        invoke(Excel.Workbooks,'Open',File);
        xlswrite1(File,v,sheet,'A1');
    catch
        %Matlab/Excel is not correctly configured
        warning('off')
        File = fullfile(dir_fig,[tit '_' sheet '.csv']);
        %Just wile a csv file
        xlswrite(File,v,sheet,'A1');
        warning('on')
    end
end
if writeET == 2 || writeET == 3
    File = fullfile(dir_fig,[tit '_' sheet '.txt']);
    fid = fopen(File,'w');
    fprintf(fid,'%1.6f  ',v);
    fclose(fid);
end