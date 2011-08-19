function data=ioi_view_time_course(IOI)

vol_anat=spm_vol(IOI.processed.file_anat);
% Basic values
nx=vol_anat(1).dim(1);
ny=vol_anat(1).dim(2);
anat_img=ioi_read_time_vol(vol_anat,1);

figure;
imagesc(anat_img');
colormap gray;
axis off;
axis equal;

h1=imrect;
wait(h1);
h2=imrect;
wait(h2);
h3=imrect;
wait(h3);


x1=round(h1.getPosition()/(2*IOI.processed.shrink_x));
x2=round(h2.getPosition()/(2*IOI.processed.shrink_x));
x3=round(h3.getPosition()/(2*IOI.processed.shrink_x));

vol_flow=spm_vol(IOI.processed.file_flow);
%vol_hbr=spm_vol(IOI.processed.file_hbr);

slice_flow=ioi_read_time_vol(vol_flow,5141);
%slice_hbr=ioi_read_time_vol(vol_hbr,5141);
% nt=length(vol_hbo);
% hbo=zeros(3,nt);
% hbr=zeros(3,nt);
% for i_slice = 1:nt
%     %We could use spm_read_vols but it brings the whole temporal file which we
%     % don't need here, instead we build our own private method
%     slice=ioi_read_time_vol(vol_hbo,i_slice);
%     tmp=slice(x1(1):x1(1)+x1(3),x1(2):x1(2)+x1(4));
%     hbo(1,i_slice)=sum(tmp(:));
%     slice=ioi_read_time_vol(vol_hbr,i_slice);
%     tmp=slice(x1(1):x1(1)+x1(3),x1(2):x1(2)+x1(4));
%     hbr(1,i_slice)=sum(tmp(:));
%     
%     tmp=slice(x2(1):x2(1)+x2(3),x2(2):x2(2)+x2(4));
%     hbo(2,i_slice)=sum(tmp(:));
%     slice=ioi_read_time_vol(vol_hbr,i_slice);
%     tmp=slice(x2(1):x2(1)+x2(3),x2(2):x2(2)+x2(4));
%     hbr(2,i_slice)=sum(tmp(:));
%     
%     tmp=slice(x3(1):x3(1)+x3(3),x3(2):x3(2)+x3(4));
%     hbo(3,i_slice)=sum(tmp(:));
%     slice=ioi_read_time_vol(vol_hbr,i_slice);
%     tmp=slice(x3(1):x3(1)+x3(3),x3(2):x3(2)+x3(4));
%     hbr(3,i_slice)=sum(tmp(:));
% 
% end
% figure;
% plot(hbo','r');
% hold on
% plot(hbr','b');
% end
data.flow=slice_flow;
%data.hbr=slice_hbr;

%%%END

