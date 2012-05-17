function IOI = extract_electro(IOI,s1,fname)
%Options:
Z.filtreFFT = 0;
Z.bits16 = 1;

sf = 10000; %sampling frequency, in Hz
fid= fopen(fname);
dir_info=dir(fname);
year=datevec(dir_info.datenum);
year=year(1); %year of the bin (used for new hardware setting in jan 2011)
if Z.bits16==1
    el=fread(fid, 'int16')'/3276.8; % conversion to 16 bits between -10 and 10
else
    el=fread(fid, 'double')';
end
%Remove stimulations (true or artefacts of stimulation)
ind = el>2; ind2 = [ind(7:end) false false false false false false];
el(ind)=el(ind2);

if s1<10,str0 = '0'; else str0=''; end
fname_el = fullfile(IOI.dir.dir_subj_res,['el_S' str0 int2str(s1) '.mat']);
save(fname_el,'el');
fclose (fid);
IOI.res.el{s1} = fname_el;
IOI.res.sfel = sf;
if Z.filtreFFT==1  %filtrage dans fourier pour enlever 60 et 5 Hz
    el=filtrageElectro33(el,sf);
end


function ve2=filtrageElectro33(ve,sf) %careful, this filter might be moving or enlarging size of 
%peaks so that electrophysiological response seems to appear to start
%before stimulation onsets

% filtrage pour acquisition à partir de 2011
% enlève le 60 Hz et le 5 Hz
f60 = 60;
f5 = 5.0031;
c1 = 0.02; %cutoff frequency band for f5 
c2 = 0.2; %cutoff frequency band for f60
cf1 = 0.3*sf; %cutoff for f5
cf2 = 0.1*sf; %cutoff for f60
f=(abs(fft(ve)));

%f=f(1:length(f));
freq=linspace(0,sf,length(f));
% plot(freq,f)
%find index of frequency to remove
index1= (((rem(freq,f5))<c1) | ((rem(freq,f5)-f5)>-c1)) & freq<cf1; %multiples of f5
index2=(((rem(freq,f60))<c2) | ((rem(freq,f60)-f60)>-c2 )) & freq<cf2; %of f60
%idem for negative fequency
index3=(((rem(freq-sf,f5))>-c1) | ((rem(freq-sf,f5)+f5)<c1)) & freq>cf1;
index4=(((rem(freq-sf,f60))>-c2) | ((rem(freq-sf,f60)+f60)<c2)) & freq>cf2;
index=index1 | index2 |index3 | index4 ;
%remove frequency component  60 and 5 Hz
f=fft(ve);
f(index)=0;
ve2=real(ifft(f));
