function [out stim couleur physio]=ioi_read_txt(fname,mode)

%vétuste
%extrait les donnees de frame.txt ou ttl.txt
%lit un ttl et extrait les valeurs importante
%path : le path du fichier .txt
% mode=1;pour le fichier frame.txt
%      =0 pour le fichier ttl.txt
% mode=2;pour le fichier info.txt
%out:[  frame temps  diffini  diff ttl couleur]
% no:le no du frame (si existant)
%sec: le temps de l'acquisition en seconde depuis le début du fichier
% diff: différence de temps entre l'acquisition et la précédente
% ttl: la valeur du ttl (si exixtant)
fid=fopen(fname);
% aa=fscanf(fid,'%c');
% aa=strread(aa,'%s','delimiter','\b');
aa=textscan(fid,'%s','delimiter','\b');
aa=aa{1};
fclose(fid);
stim=   zeros(length(aa),1);
frame=   zeros(length(aa),1);
couleur=   char(zeros(length(aa),1));
sec=   zeros(length(aa),1);
ttl=   zeros(length(aa),1);
diff=   zeros(length(aa),1);
tini=0;
physio.Resp=zeros(length(aa),1);
physio.bpm=zeros(length(aa),1);
physio.CO2=zeros(length(aa),1);
physio.Temperature=zeros(length(aa),1);
ajoutframe=0;
z=0;
for i1=1:length(aa)
    y=i1;
    tline = aa{y};
    if tline==-1,
        Disp2(handles,['readTXT: Fichier vide : ' path]);
        break;
    end
    if mode==1
        frame1=sscanf(tline(8:13),'%f');
        if frame1(end)<0 && framevieux(end)>0;
            ajoutframe=ajoutframe+65536;
        end
        frame(y,1)= frame1(end)+ajoutframe;
        framevieux=frame1;
        matches = findstr(tline, 'Temps :');
        matches=matches+7;
    elseif mode==2
        tline(40)=' ';
        if ~strcmp(tline(1:4),'Bloc')
            continue
        end
        z=z+1;
        y=z;
        Bloc1(y)=sscanf(tline(7:13),'%f');
        matches = findstr(tline, 'Frame ini :');
        frame(y)=sscanf(tline(matches+11:matches+15),'%f');
        matches = findstr(tline, 'Temps :');
        matches=matches+7;
    elseif mode==0
        frame(y,1)=-1;
        matches=1;
    end
    match2=findstr(tline, 'Couleur :');
    if ~isempty(match2)&&length(tline)>=match2+10
        couleur(y,1)=tline(match2+10);
    else
        couleur(y,1)='-';
    end
    match2=findstr(tline, 'Resp :');
    if ~isempty(match2)&&length(tline)>=match2+9
        physio.Resp(y,1)=sscanf(strrep(tline(match2+7:match2+9),',','.'),'%f');
    else
        physio.Resp(y,1)=NaN;
    end
    match2=findstr(tline, 'bpm :');
    if ~isempty(match2)&&length(tline)>=match2+12
        physio.bpm(y,1)=sscanf(strrep(tline(match2+6:match2+12),',','.'),'%f');
    else
        physio.bpm(y,1)=NaN;
    end
    match2=findstr(tline, 'Press :');
    if ~isempty(match2)&&length(tline)>=match2+12
        aaa=strrep(tline(match2+7:match2+12),',','.');
        physio.Press(y,1)=sscanf(strrep(tline(match2+7:match2+12),',','.'),'%f');
    else
        physio.Press(y,1)=NaN;
    end
    match2=findstr(tline, 'Tempe :');
    if ~isempty(match2)&&length(tline)>=match2+11
        physio.Temperature(y,1)=sscanf(strrep(tline(match2+7:match2+11),',','.'),'%f');
    else
        physio.Temperature(y,1)=NaN;
    end
    match2=findstr(tline, 'CO2 :');
    if ~isempty(match2)&&length(tline)>=match2+10
        physio.CO2(y,1)=sscanf(strrep(tline(match2+6:match2+10),',','.'),'%f');
    else
        physio.CO2(y,1)=NaN;
    end
    %calcul différente valeur de temps
    time2=tline(matches+11:matches+18);
    milisec=sscanf(tline(matches+20:matches+22),'%f')/1000;
    time2=[ 0 0 0 sscanf(time2,'%f:')'];
    vecteur = [0 0 0 time2(4:end) milisec] ;
    sec(y,1)=vecteur(4)*3600+vecteur(5)*60+vecteur(6)+vecteur(7);
    if y==1
        tini=vecteur(4)*3600+vecteur(5)*60+vecteur(6)+vecteur(7);
    else
        diff(y,1)=sec(y,1)-sec(y-1,1);
    end
    %correctionseuil(aa,vecteur,y,diff) % à laisser peut encore être utile
    if mode==0
        finDeLigne=tline(matches+23:end);
        if strcmp(finDeLigne(1:4),' , S'), finDeLigne=finDeLigne(10:end); end
        finDeLigne=sscanf(regexprep(finDeLigne, ',', '.'),'%f');
        ttl(y,1)=finDeLigne;
    elseif mode==2
        ttl(y,1)=Bloc1(y);
    elseif mode==1
        ttl(y,1)=0;
    end
    %     % A CONSERVER, utile pour appeler la fonction corrigeframetxt.m
    %     if y==1, fid=fopen([path(1:end-4) '2.txt'],'w');end
    %     corrigeframetxt([path(1:end-10) '_images'],'LL',frame(y,1),tline,fid);
    %     if y== length(aa),  fclose(fid);    end
end
physio.CO2=single(physio.CO2);
physio.bpm=single(physio.bpm);
physio.Resp=single(physio.Resp);
physio.Press=single(physio.Press);
physio.frame=single(frame);
physio.temps=single(sec-tini);
physio.Temperature=single(physio.Temperature);
out=[frame sec sec-tini diff ttl];
out=out(out(:,1)~=0,:);
if nargout==2 && mode==0 %détermine à quelle stimulation on est rendu
    index0=find(out(:,4)>0);
    index=find(out(:,4)>2.1);
    if ~isempty(index0) && index0(1)~=index(1) % si la première stimulation arrive au temps 0, permet de la prendre en compte
        index=[index0(1); index];
    end
    if isempty(min(index)>10) % si on a pas de stim, permet d'analyser quand même
        index=1;
    end
    if size(out,1)>=2
         stim =[out(index,1:4) round(out(min(index+1,size(out,1)),5)*50)/500] ; %arrondi et met en mV
    else
        stim =[out(index,1:4) 0] ;
    end
end
