function IOI = disp_msg(IOI,msg)
disp(msg);
try
    IOI.warning = [IOI.warning; msg];
catch
    IOI.warning = [IOI.warning msg];
end