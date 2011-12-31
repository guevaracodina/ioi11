function frame = ioi_check_frame(frame,nF)
frame = round(frame);
if frame > nF
    frame = nF;
end
if frame < 1;
    frame = 1;
end
