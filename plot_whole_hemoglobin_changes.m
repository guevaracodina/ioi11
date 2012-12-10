clc
clear all
load('W:\CONG_TRI\Res\12_11_19,CO22\ROI.mat');
load('W:\CONG_TRI\Res\12_11_19,CO22\IOI.mat');
Gaussion_LP_filter=0;
LPF.fwhm1=0.5;
data=ROI;
IOI.dev.TR=0.2;
fs=5;
lpf_hemoglobin_butter_On=1;
HPF.hpf_butter_On=1;
HPF.hpf_butter_freq=0.1;
HPF.hpf_butter_order=3;
lpf_butter_On=0;
cutoff=0.2;
%         cutoff=0.025;
FilterOrder=3;
for s=12
    for r=1:length(ROI)
        if HPF.hpf_butter_On
            for i=5:7
            tmp_DC = median(data{1,r}{s,i});
            data{1,r}{s,i} = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,data{1,r}{s,i});
            data{1,r}{s,i} = data{1,r}{s,i} + tmp_DC; %add back the DC component
            end
        end

%         if Gaussion_LP_filter
%             for i=5:7
%                 K = get_K(1:length(data{1,r}{s,i}),LPF.fwhm1,IOI.dev.TR);
%                 data{1,r}{s,i} = ioi_filter_HPF_LPF_WMDL(K,(data{1,r}{s,i})')';
%             end
%         end
        
        if lpf_hemoglobin_butter_On
            for i=5:7
                data{1,r}{s,i}=ButterLPF(fs,cutoff,FilterOrder,data{1,r}{s,i});
            end
        end
        
        
        HbO{r}=data{1,r}{s,5};
        HbO{r}=HbO{r}-mean(HbO{r})+60;
        HbO{r}=100*(HbO{r}/60-1);
%         mean_HbO=mean(HbO{r});
%         HbO{r}=(HbO{r}-mean_HbO)./mean_HbO;
        a=figure;
        lp=linspace(0,length(HbO{r})/fs,length(HbO{r}));
        plot(lp,HbO{r},'r'); hold on;
        
        HbR{r}=data{1,r}{s,6};
        HbR{r}=HbR{r}-mean(HbR{r})+40;
        HbR{r}=100*(HbR{r}/40-1);
%         mean_HbR=mean(HbR{r});
%         HbR{r}=(HbR{r}-mean_HbR)./mean_HbR;
%         lp=linspace(0,length(HbR{r})/fs,length(HbR));
        plot(lp,HbR{r},'b'); hold on;
        
        HbT{r}=data{1,r}{s,5}+data{1,r}{s,6};
        HbT{r}=HbT{r}-mean(HbT{r})+100;
        HbT{r}=100*(HbT{r}/100-1);
%            mean_HbT=mean(HbT{r});
%         HbT{r}=(HbT{r}-mean_HbT)./mean_HbT;
        
        Flow{r}=data{1,r}{s,7};
%         cutoff_flow=0.08;
%         FilterOrder_flow=3;
%         if lpf_butter_On            
%             Flow=ButterLPF(fs,cutoff_flow,FilterOrder_flow,Flow);
% %         end
%         %********CO14 session 05
%         x_line=[1 154 365 518 776];
%         x1=[13 183 395 553 810];
%         %**********CO14 session 06
%         x_line=[400];
%         x1=[428];
%*********************end
%**********************CO22 session7
% x_line=[317 676];
% x1=[400 760];
%***************end
% %**********cO22 session8
% x_line=[209 548];
% x1=[290 621];
%********end

% %**********cO22 session9
% x_line=[40 302 565 850];
% x1=[132 384 655 1058];
% %********end
% %**********cO22 session10
% x_line=[254 504 880];
% x1=[315 575 940];
% %********end
%**********cO22 session12
x_line=[65 308];
x1=[121 515];
%********end

%         lpf_butter_On=0;
%         cutoff=0.05;
%         Flow = ButterLPF(fs,cutoff,FilterOrder,Flow);
%         x=[8 132 362 481 580 726]; x1=[47 166 400 514 614 762];
% %         x=[0 155 367 520 780];x1=[13 183 395 553 810];
        b=mean(Flow{r});
        Flow{r}=100*((Flow{r}-b)./b);
        plot(lp,HbT{r},'g');hold on;
        plot(lp,Flow{r},'k');hold on;
%         x=[HbO' HbR' HbT' Flow'];
%         ORTF_correlation=[];
%         ORTF_correlation=corr(x);
        
%         lp=linspace(0,length(HbT)/fs,length(HbT));
        
        for i=1:length(x_line)
        plot([x_line(i),x_line(i)],[min(HbR{r}),max(HbO{r})],'m');hold on;
        plot([x1(i),x1(i)],[min(HbR{r}),max(HbO{r})],'c');hold on;
        end
%         [AX, H1,H2]=plotyy(lp,HbT,lp,Flow); hold on;
%         set(H1,'Color','g'); set(H2,'Color','k');
% %         plot(lp,Flow,'k'); hold on;
%       set(get(AX(1),'Ylabel'),'String','percent changes of HbO, HbR and HbT');
%       set(get(AX(2),'Ylabel'),'String','percent changes of Flow');
        
        xlabel('time(s)');
        ylabel('percent changes');
%         ylabel('percent change');
        tit = ['Session ' int2str(s) ', ROI ' IOI.ROIname{r}];
        title(tit);
%         legend('HbO','HbR','HbT','onsets','the end of seizures');
        legend('HbO','HbR','HbT','Flow','seizure onset','seizure offset');
        filename=['W:\CONG_TRI\Res\12_11_19,CO22\each seizure\s12\',tit];
        saveas(a,filename, 'jpg');
        saveas(a,filename, 'fig');
        close(a);        
    end
end
hbo=HbO;
hbt=HbT;
hbr=HbR;
flow=Flow;
s=5;
% %************CO14 session 5
x_line=[1 130 340 500 750]; %onset-10s
x1=[13 183 395 553 810];
% %********end
        %**********CO14 session 06
%         x_line=[390];
%         x1=[428];
% % *********************end
%***************plot all hemoglobin changes in all ROI during one seizure 
for i=1:length(x_line)
    %**************plot HbO
e=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for r=1:length(ROI)
   lp=linspace(0,x1(i)-x_line(i),length(HbO{r}(x_line(i)*fs:x1(i)*fs)));
   plot(lp,HbO{r}(x_line(i)*fs:x1(i)*fs));   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['HbO changes for all areas in', ' Session ' int2str(s) ,' seizure' int2str(i)];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all areas in one seizure\',tit];
saveas(e,filename, 'jpg');
saveas(e,filename, 'fig');
close(e);

%*************plot HbR
b=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for r=1:length(ROI)
   lp=linspace(0,x1(i)-x_line(i),length(HbR{r}(x_line(i)*fs:x1(i)*fs)));
   plot(lp,HbR{r}(x_line(i)*fs:x1(i)*fs));   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['HbR changes for all areas in', ' Session ' int2str(s) ,' seizure' int2str(i)];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all areas in one seizure\',tit];
saveas(b,filename, 'jpg');
saveas(b,filename, 'fig');
close(b);
%***************plotHbT
c=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for r=1:length(ROI)
   lp=linspace(0,x1(i)-x_line(i),length(HbT{r}(x_line(i)*fs:x1(i)*fs)));
   plot(lp,HbT{r}(x_line(i)*fs:x1(i)*fs));   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['HbT changes for all areas in', ' Session ' int2str(s) ,' seizure' int2str(i)];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all areas in one seizure\',tit];
saveas(c,filename, 'jpg');
saveas(c,filename, 'fig');
close(c);

%******plot flow changes

d=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for r=1:length(ROI)
   lp=linspace(0,x1(i)-x_line(i),length(Flow{r}(x_line(i)*fs:x1(i)*fs)));
   plot(lp,Flow{r}(x_line(i)*fs:x1(i)*fs));   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['Flow changes for all areas in', ' Session ' int2str(s) ,' seizure' int2str(i)];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all areas in one seizure\',tit];
saveas(d,filename, 'jpg');
saveas(d,filename, 'fig');
close(d);

end

%***************plot all seizure changes in one ROI 
%can not do it because the duration was different

HbO=hbo;
HbR=hbr;
HbT=hbt;
Flow=flow;
% s=5;
% % %************CO14 session 5
% x_line=[1 154 365 518 776]; %onset-10s
% x1=[13 183 395 553 810];
% % %********end
        %**********CO14 session 06
%         x_line=[390];
%         x1=[428];
% *********************end
%***************plot all hemoglobin changes in all ROI during one seizure 
for i=1:length(x_line)
    duration(i)=x1(i)-x_line(i);
end
duration_mean=fix(mean(duration));
for r=1:length(ROI)
    %**************plot HbO
e=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for j=1:length(x_line)
    HbO_seizure{j}=resample(HbO{r}(x_line(j)*fs:x1(j)*fs), duration_mean*fs,length(HbO{r}(x_line(j)*fs:x1(j)*fs))); 
   lp=linspace(0,duration_mean,length(HbO_seizure{j}));
   plot(lp,HbO_seizure{j});   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['HbO changes for all seizures in', ' Session ' int2str(s) ,' ROI' IOI.ROIname{r}];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all seizures in one area\',tit];
saveas(e,filename, 'jpg');
saveas(e,filename, 'fig');
close(e);

%*************plot HbR
b=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for j=1:length(x_line)
    HbR_seizure{j}=resample(HbR{r}(x_line(j)*fs:x1(j)*fs), duration_mean*fs,length(HbR{r}(x_line(j)*fs:x1(j)*fs))); 
   lp=linspace(0,duration_mean,length(HbR_seizure{j}));
   plot(lp,HbR_seizure{j});   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['HbR changes for all seizures in', ' Session ' int2str(s) ,' ROI' IOI.ROIname{r}];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all seizures in one area\',tit];
saveas(b,filename, 'jpg');
saveas(b,filename, 'fig');
close(b);
%***************plotHbT
c=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for j=1:length(x_line)
    HbT_seizure{j}=resample(HbT{r}(x_line(j)*fs:x1(j)*fs), duration_mean*fs,length(HbT{r}(x_line(j)*fs:x1(j)*fs))); 
   lp=linspace(0,duration_mean,length(HbT_seizure{j}));
   plot(lp,HbT_seizure{j});   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['HbT changes for all seizures in', ' Session ' int2str(s) ,' ROI' IOI.ROIname{r}];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all seizures in one area\',tit];
saveas(c,filename, 'jpg');
saveas(c,filename, 'fig');
close(c);

%******plot flow changes

d=figure;
ColorSet = varycolor(length(ROI));
set(gca, 'ColorOrder', ColorSet);
hold all
for j=1:length(x_line)
    Flow_seizure{j}=resample(Flow{r}(x_line(j)*fs:x1(j)*fs), duration_mean*fs,length(Flow{r}(x_line(j)*fs:x1(j)*fs))); 
   lp=linspace(0,duration_mean,length(Flow_seizure{j}));
   plot(lp,Flow_seizure{j});   
end
hc1 = ioi_set_colorbar(gcf,length(ROI));
legend off
set(gcf, 'Colormap', ColorSet);
tit = ['Flow changes for all seizures in', ' Session ' int2str(s) ,' ROI' IOI.ROIname{r}];
title(tit);
filename=['W:\CONG_TRI\Res\12_08_27,CO14_new\ROI_for_fred\each seizure\s05\all seizures in one area\',tit];
saveas(d,filename, 'jpg');
saveas(d,filename, 'fig');
close(d);

end





%************************ plot HbO HbR and flow

% clear all
% 
% try
%     %clear all;
%     clc;
%     %************************ plot HbO HbR and flow
% load('W:\CONG_TRI\Res\12_08_27,CO14_old\ROI_seizure\ROI.mat');
% load('W:\CONG_TRI\Res\12_08_27,CO14_old\ROI_seizure\IOI.mat');
%     %     fs
%     %      Y = ButterLPF(fs,cutoff,FilterOrder,Y);
%     %      LPF = ioi_get_LPF(job);
%     % tmp_pars = IOI.sess_res{5}.parameters;
%     
%     fs=5;
% %      tb{1}=8;ta{1}=47;
% %     tb{2}=132;ta{2}=166;
% %     tb{3}=362;ta{3}=400;
% %     tb{4}=481;ta{4}=514;
% %     tb{5}=580;ta{5}=614;
% %     tb{6}=726;ta{6}=762;
% %      tb{1}=8;ta{1}=114;
% %     tb{2}=257;ta{2}=289;
% %     tb{3}=344;ta{3}=431;
% %     tb{4}=600;ta{4}=724;
% % tb{1}=1;ta{1}=13;
% tb{1}=155;ta{1}=183;
% tb{2}=367;ta{2}=395;
% tb{3}=520;ta{3}=553;
% tb{4}=780;ta{4}=810;
%     Nz=4;
%     %c0=5;    % co is the number of session
%     c0=5;
%     time_baseline=10;
%     time_afterSZ=10;
%     time_offset=5;
%     resampleOn = 1; %if don't resample,  truncate to shortest seizure
%     show_percent_changes = 1;
%     %     cutoff=130;
%     %     FilterOrder=3;
%     %     data{1,r1}{c0,j} = ButterLPF(fs,cutoff,FilterOrder,data{1,r1}{c0,j});
%     LPF.fwhm1=0.5;
%     IOI.dev.TR=0.2;
%     Gaussion_LP_filter=1;
%     HG=cell(1,4);
% %     for i=1:4
% %         HG{i}=struct();
% %     end
%     HG{1}.name1='HbO';
%     HG{2}.name1='HbR';
%     HG{3}.name1='Flow';
%     HG{4}.name1='HbT';
%     HG{1}.name2='HbO_baseline';
%     HG{2}.name2='HbR_baseline';
%     HG{3}.name2='Flow_baseline';
%     HG{4}.name2='HbT_baseline';
%     HG{1}.name3='HbO_seizure';
%     HG{2}.name3='HbR_seizure';
%     HG{3}.name3='Flow_seizure';
%     HG{4}.name3='HbT_seizure';
%     HG{1}.name2_mean='HbO_mean';
%     HG{2}.name2_mean='HbR_mean';
%     HG{3}.name2_mean='Flow_mean';
%     HG{4}.name2_mean='HbT_mean';
%     
%     data=ROI;
%     for r1=1:length(ROI)
%         %***************filter the data
%         if Gaussion_LP_filter
%             for i=5:7
%                 K = get_K(1:length(data{1,r1}{c0,i}),LPF.fwhm1,IOI.dev.TR);
%                 data{1,r1}{c0,i} = ioi_filter_HPF_LPF_WMDL(K,(data{1,r1}{c0,i})')';
%             end
%         end
%         %**caculate the baseline
%         for i=1:Nz   %%%%%%%%%%%%%the number of seizure
%             index1{i}=(tb{i}*fs-time_baseline*fs-time_offset*fs):(tb{i}*fs-time_offset*fs);
%             index{i}=(tb{i}*fs-time_offset*fs):(ta{i}*fs+time_afterSZ*fs);
%             sz_duration(i)=length(index{i});
%             for j=1:4          %**************HbO HbR Flow and HbT
%                 switch j
%                     case 1 %HbO
%                         HG{j}.value1=data{1,r1}{c0,5};
%                         HG{j}.value1=HG{j}.value1-mean(HG{j}.value1)+60; %normalization
%                         if show_percent_changes
%                             HG{j}.value1=100*(HG{j}.value1/60-1);
%                         end
%                         HG{j}.value2(i,:)= HG{j}.value1(index1{i});    %%%%%%%%%%%%%calculate the baseline
%                         HG{j}.value3{i}= HG{j}.value1(index{i});     %%%%%%%%%%%%%%%%%%%calculate the seizure
%                     case 2%*HbR
%                         HG{j}.value1=data{1,r1}{c0,6};
%                         HG{j}.value1=HG{j}.value1-mean(HG{j}.value1)+40;
%                         if show_percent_changes
%                             HG{j}.value1=100*(HG{j}.value1/40-1);
%                         end
%                         HG{j}.value2(i,:)= HG{j}.value1(index1{i});
%                         HG{j}.value3{i}= HG{j}.value1(index{i});
%                     case 3  %flow
%                         HG{j}.value1=data{1,r1}{c0,7};
%                         HG{j}.value2(i,:)= HG{j}.value1(index1{i});
%                         HG{j}.value3{i}= HG{j}.value1(index{i});
%                     case 4   %HbT
%                         HG{j}.value1=(data{1,r1}{c0,5}+data{1,r1}{c0,6})-mean((data{1,r1}{c0,5}+data{1,r1}{c0,6}))+100;
%                         if show_percent_changes
%                             HG{j}.value1=100*(HG{j}.value1/100-1);
%                         end
%                         HG{j}.value2(i,:)= HG{j}.value1(index1{i});
%                         HG{j}.value3{i}= HG{j}.value1(index{i});
%                 end
%                 HG{j}.value2_mean=mean(HG{j}.value2,2);  %%%%%%%%%%%%%%%%calculation the average of baseline
% 
%             end
%         end
%             for i=1:Nz
%                 b=figure;
%                     %                     HG{j}.value3{i}= resample(HG{j}.value3{i},duration_mean,length(HG{j}.value3{i}));
%                     HG{1}.value4(i,:)=HG{1}.value3{i}-HG{1}.value2_mean(i);
%                     HG{2}.value4(i,:)=HG{2}.value3{i}-HG{2}.value2_mean(i);
%                     HG{4}.value4(i,:)=HG{4}.value3{i}-HG{4}.value2_mean(i);
%                     lp=linspace(-time_offset,length(HG{1}.value4)/fs-time_offset,length(HG{1}.value4));
%                     plot(lp,HG{1}.value4,'r'); hold on;
%                     plot(lp,HG{2}.value4,'b'); hold on;
%                     plot(lp,HG{4}.value4,'g'); hold on;
%                     xlim([-time_offset length(HG{1}.value4)/fs-time_offset]);
%                     xlabel('time(s)');
%                     ylabel('percent change');
%                     tit = ['Session ' int2str(c0) ', ROI ' IOI.ROIname{r1}  ', seizure ' int2str(i)] ;
%                     title(tit);
%                     legend('HbO','HbR','HbT');
%                     filename=['W:\CONG_TRI\Res\12_08_27,CO14_old\ROI_seizure\each seizure\s05\',tit];
%                     saveas(b,filename, 'tif');
%                     saveas(b,filename, 'fig');
%                     close (b);
%                     HG{1}.value4=[];
%                     HG{2}.value4=[];
%                     HG{4}.value4=[];
%             end
% 
%     end
% catch exception
%     disp(exception.identifier)
%     disp(exception.stack(1))
% end
% 
