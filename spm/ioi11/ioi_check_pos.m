function [xpos ypos] = ioi_check_pos(xpos,ypos,nx,ny)
xpos = round(xpos);
ypos = round(ypos);
if xpos > nx
    xpos = nx;
end
if xpos < 1
    xpos = 1;
end
if ypos > ny
    ypos = ny;
end
if ypos < 1
    ypos = 1;
end

end