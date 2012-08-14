function [sts B] = ioi_check_color_order(image_total,str1,str_laser)
switch str1
    case 'R'        
        thr = 10; %threshold
    case 'G'
        thr = 10; %threshold
    case 'Y'
        thr = 10; %threshold
    case 'L'
        thr = 1; %threshold
end
[nx ny nn nt] = size(image_total);
nx1 = round(nx/3);
nx2 = round(2*nx/3);
ny1 = round(ny/2);
%ny2 = round(2*ny/3);
bfr = [];
bt = [];
%criterion
sts = 1;
for i0 = 1:nt
    t = std(image_total(nx1:nx2,ny1,:,i0)-image_total((nx1-1):(nx2-1),ny1,:,i0))/mean(abs(image_total(nx1:nx2,ny1,:,i0)));
    if strcmp(str1,str_laser)
        if t < thr %not enough variance for laser speckle
            sts = 0;
            bfr = [bfr i0];
            bt = [bt t];
        end
    else
        if t > thr %too much variance for green, yellow or red colors
            sts = 0;
            bfr = [bfr i0];
            bt = [bt t];
        end
    end
end
B.bfr = bfr;
B.bt = bt;

% while i0 < nt
%     i0 = i0+1;
%     t = std(image_total(nx1:nx2,ny1,:,i0)-image_total((nx1-1):(nx2-1),ny1,:,i0))/mean(image_total(nx1:nx2,ny1,:,i0));
%     if strcmp(str1,str_laser)
%         if t < thr %not enough variance for laser speckly
%             sts = 0;
%             break
%         end
%     else
%         if t > thr %too much variance for green, yellow or red colors
%             sts = 0;
%             break
%         end
%     end
% end
