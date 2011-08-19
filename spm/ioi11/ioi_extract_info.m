function [info physio]=ioi_extract_info(file_resume,file_info,file_TTL,file_Frame)
physio=[];  
info=[];
% Extract experiment resume
%file_resume = [files_directory, filesep, '..' filesep, 'experiment_resume.txt'];
info=private_doResume(file_resume,info);

% Extract scan info
%file_info = [files_directory, filesep, 'info.txt'];
info=private_doInfo(file_info,info);

% Extract Stimulation info
[info physio]=private_doTTL(file_info,file_TTL,file_Frame, info);

%faire physio si nécessaire
info=private_doPhysio(info,physio);


% Read file 'resume_manip.txt'
% Internal function, does not need to be accessed elsewhere
function info=private_doResume(fname,info)

info.age=0;
info.sexe=[];
info.ne=[];
info.numero='AA00E00';
if exist(fname,'file')
    fid=fopen(fname);
    aa=textscan(fid,'%s','delimiter','\b','whitespace', '');
    aa=aa{1};
    fclose(fid);
    % Initializations
    info.numero=' ';
    info.poids=' ';
    info.seuil10ms=' ';
    info.seuil1ms=' ';
    info.date=' ';
    info.rate=' ';
    info.rein=' ';
    info.foie=' ';
    info.seuiltreshold=' '; 
    info.ne= ' ';  
    info.age= ' ';
    % Now reading the file
    for y=1:length(aa)
        tline = aa{y};
        k=strfind( tline,':');
        if strcmpi(tline(1:min(8,length(tline))),'numéro :')
            info.numero=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(6,length(tline))),'poid :')
            info.poids=strtok(tline(k+1:end));
        elseif any(strcmpi(tline(1:min(4,length(tline))),{'Né :','ne :','Ne :','né :'}))
            info.ne=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(7,length(tline))),'poids :')
            info.poids=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(15,length(tline))),'seuil à 10 ms :')
            info.seuil10ms=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(15,length(tline))),'seuil à  1 ms :') || strcmpi(tline(1:min(13,length(tline))),'seuil à 1 ms :')
            info.seuil1ms=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(6,length(tline))),'date :')
            info.date=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(6,length(tline))),'rate :')
            info.rate=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(6,length(tline))),'rein :')
            info.rein=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(6,length(tline))),'sexe :')
            info.sexe=strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(6,length(tline))),'foie :')
            info.foie =strtok(tline(k+1:end));
        elseif strcmpi(tline(1:min(16,length(tline))),'seuil treshold :')
            info.seuiltreshold=strtok(tline(k+1:end));
        end
    end
    % Find age, conversion in months
    if length(info.ne)>3
        info.age=datevec(datenum(info.date, 'yy_mm_dd', 1400)-datenum(info.ne, 'yy_mm_dd', 1400));
        info.age= [info.age(1)*12 + (info.age(2)) + (info.age(3))/30];
    end
    
    % Experiment number, priorities as follows: read it from the file, if
    % not present, take it from the filename, otherwize set to 'z'
    numero1=fname(end+(-7:-1));
    [numero1 exp]=strtok(numero1,'E');
    if isfield(info,'numero')
        numero=info.numero;
    else
        if ~isempty(numero1)
            numero=numero1;
        else
            numero='z';
        end
    end
    if ~isempty(exp)
        info.numero=[numero exp];
    end
    
else
    disp('private_doResume: Could not find resume_manip.txt')
end

% Read file info.txt, does not need to be used elsewhere
function info=private_doInfo(fname,info)

fid=fopen(fname);
aa=textscan(fid,'%s','delimiter','\b','whitespace', '');
aa=aa{1};
fclose(fid);

y = 1;
y2 = 1;
info.commentaire='';
info.stim=0;
info.bloc=[];
info.Y=0;
info.dY=0;
info.X=0;
info.dX=0;
info.dtexposition=0;
info.dtframe=0;
info.lumiere=' ';
info.protocolelumiere=0;
info.commentaire=' ';
info.protocoleStimulation =' ';
info.tempsIni='';
info.olddY=[];
info.olddX=[];
info.stimulateur='';
for y=1:length(aa)
    tline = aa{y};
    k=strfind( tline,' : ');
    if strcmpi(tline(1:min(1,length(tline))),'Y')
        info.Y=sscanf(strrep(tline(5:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(9,length(tline))),'Hauteur Y')
        info.dY=sscanf(strrep(tline(12:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(12,length(tline))),'oldHauteur Y')
        info.olddY=sscanf(strrep(tline(15:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(12,length(tline))),'oldLargeur X')
        info.olddX=sscanf(strrep(tline(15:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(1,length(tline))),'X')
        info.X=sscanf(strrep(tline(5:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(9,length(tline))),'Largeur X')
        info.dX=sscanf(strrep(tline(12:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(13,length(tline))),'Exposition ms')
        info.dtexposition=sscanf(strrep(tline(16:end),',','.'),'%f')/1000;
    elseif strcmpi(tline(1:min(23,length(tline))),'periode sauvegarde (ms)')
        info.dtframe=sscanf(strrep(tline(26:end),',','.'),'%f')/1000;
    elseif strcmpi(tline(1:min(7,length(tline))),'Lumière')
        info.lumiere=tline(10:end);
        info.lumiere=strrep(info.lumiere,'Laser','L');
    elseif strcmpi(tline(1:min(17,length(tline))),'Protocole lumière')
        info.protocolelumiere=tline(20:end);
    elseif strcmpi(tline(1:min(11,length(tline))),'Commentaire')
        ligneCommentaire=y;
        commentaire={tline(14:end)};
    elseif strcmpi(tline(1:min(22,length(tline))),'Protocole stimulation ')
        info.protocoleStimulation =tline(25:end);
    elseif strcmpi(tline(1:min(14,length(tline))),'Seuil Stim, mA')
        info.stim =sscanf(strrep(tline(17:end),',','.'),'%f');
    elseif strcmpi(tline(1:min(11,length(tline))),'Stimulateur')
        info.stimulateur =tline(14:end);
    elseif strcmpi(tline(1:min(7,length(tline))),'Bloc :0')
        if isempty(info.tempsIni)
            match=strfind(tline,'Temps');
            info.tempsIni=tline(match+18:match+25);
        end
    end
    if    strcmp( tline(1:min(4,length(tline)))  ,'Path')
        info.Path=tline(8:end);
        linePath=y;
    end
    if ~exist('linePath','var') && exist('ligneCommentaire','var') && y~=ligneCommentaire  % si la signe de commentaire est sur plusieurs lignes
        commentaire=[commentaire ; {tline}];
    end
    %  pour les bloc de stimulation et leur temps respectifs
    if strcmp( tline(1:min(4,length(tline)))  ,'Bloc')
        frame(y2,1)=sscanf(tline(7:14),'%f');
        matches = findstr(tline, 'Temps :');
        matches=matches+7;
        date(1)=str2num(tline(matches:matches+3));
        date(2)=str2num(tline(matches+5:matches+6));
        date(3)=str2num(tline(matches+8:matches+9));
        time=tline(matches+11:matches+18);
        milisec=str2num(tline(matches+20:matches+22))/1000;
        finDeLigne=tline(matches+23:end);
        finDeLigne=sscanf(regexprep(finDeLigne, ',', '.'),'%f');
        time=datenum(time, 'HH:MM:SS');
        time=datevec(time);
        vecteur=[date(1:3) time(4:end) milisec];
        sec(y2,1)=vecteur(4)*3600+vecteur(5)*60+vecteur(6)+milisec;
        y2=y2+1;
    end
end
info.commentaire=char(commentaire);
if exist('vecteur','var')
    temps=vecteur;
else
    temps=[];
end
if (~isfield(info,'stim')|| info.stim==0 || info.stim>.4) && (isempty(info.stimulateur) || strncmp(info.stimulateur,'5 V',3) )
    info.stimvieux=info.stim;
    text=info.commentaire';
    for i1=1:length(text); %trouver l'intensité de stimulation dans le fichier info.txt
        stim=sscanf(text(i1:end),'%f');
        if isnumeric(stim)
            info.stim=stim;
            break
        end
    end
end

 
function [info physio]=private_doTTL(file_info,file_TTL,file_Frame,info)

% boucle pour lire les infos du TTL.txt
info.stim1=[]; 
info.stim2=[];
info.Tacq=0;

%extrait les info sur ls frames
[info.Frame info.FrameCouleur info.FrameCouleur physio]=ioi_read_txt(file_Frame,1); 
 
if size(info.Frame,1)~=0
    [info.bloc]=ioi_read_txt(file_info,2); %
    if isempty(info.bloc), info.bloc=[ 0 0 0 0 0]; end
    [info.ttl info.stim1 ]=ioi_read_txt(file_TTL,0);
    % corection pour alligner exactement les stim aux Frame (11 aout 09) corrige d'environ .1 s
    info.stim1(:,3)=info.stim1(:,2)-info.Frame(1,2);
    if size(info.stim,1)>1 &&( (info.bloc(end,2)-info.stim1(end-1,2))>19 && (info.stim1(end,2)-info.bloc(end,2))>3  )
        info.stim1(end,2:4)=info.bloc(end,2:4);
    end
    [info.Frame info.stim2]=framestim(info.Frame, info.stim1);
    info.Tacq=info.Frame(end,3)/info.Frame(end,1);
    %déterminer la stim arrive dans quel bloc, et à quelle position dans le
    %bloc
    for i1=1:size(info.stim1,1)
        info.stim1(i1,6)=find(info.stim1(i1,3)>(info.bloc(:,3)-3),1,'last');
        info.stim1(i1,7)=sum(info.stim1(1:i1,6)==info.stim1(i1,6) & info.stim1(1:i1,5)~=0);
    end
    if length(unique(info.stim1(:,6)))==1 %corriger un fuck des premiers fichiers ou on utilisait pas la structure par bloc
        info.stim1=info.stim1(:,[ 1 2 3 4 5 7 6]);
    end
    %déterminer les valeurs des différentes stimulations
    info.stim3=unique(info.stim1(:,7));
    info.stim3=info.stim3(info.stim3~=0);
    for i2=1:length(info.stim1)
        index=find(info.ttl(:,2)==info.stim1(i2,2));
        info.ttl(index(1):end,6)=info.stim1(i2,7);
    end
    for i1=1:size(info.stim3,1)
        value=sort(info.stim1(info.stim1(:,7)==info.stim3(i1),5));
        info.stim3(i1,2)=median(value); %intensité mediane
        aa=info.ttl(info.ttl(:,6)==i1,:);
        aa=aa(aa(:,4)<2 & aa(:,4)~=0,:);
        info.stim3(i1,3)=mean(aa(:,4)); %intervale entre 2 stim électrique
        info.stim3(i1,4)=round(size(aa,1)/sum(info.stim1(:,7)==i1))+1; %nbre stim électrique
        info.stim3(i1,5)= info.stim3(i1,3) *info.stim3(i1,4); %durée totale d'un bloc de stim
    end
    physio.stim2=info.stim2; 
end




function [info]=private_doPhysio(info,physio)
% détermine les statistiques du fichier physio
if exist('physio','var')
    info.physio.resp=mean(physio.Resp);    
    info.physio.respstd=std(physio.Resp);
    info.physio.press=mean(physio.Press);
    info.physio.pressstd=std(physio.Press);
    info.physio.CO2=mean(physio.CO2);
    info.physio.CO2std=std(physio.CO2);
    aa=reshape(physio.CO2(1:floor(length(physio.CO2)/20)*20),20,floor(length(physio.CO2)/20));
    info.physio.CO2exp=mean(max(aa,[],1));
    info.physio.CO2expstd=std(max(aa,[],1));
    info.physio.bpm=mean(physio.bpm);
    info.physio.bpmstd=std(physio.bpm);
    info.physio.Temperature=mean(physio.Temperature);
    info.physio.Temperaturestd=std(physio.Temperature);
    info.bpm=[ num2str(info.physio.bpm,'%4.1f') ' +- ' num2str(info.physio.bpmstd,'%4.1f') ];
    info.press=[ num2str(info.physio.press,'%4.1f') ' +- ' num2str(info.physio.pressstd,'%4.1f') ];
    info.CO2=[ num2str(info.physio.CO2,'%4.2f') ' +- ' num2str(info.physio.CO2std,'%4.2f') ];
    info.CO2exp=[ num2str(info.physio.CO2exp,'%4.2f') ' +- ' num2str(info.physio.CO2expstd,'%4.2f') ];
    info.resp=[ num2str(info.physio.resp,'%4.3f') ' +- ' num2str(info.physio.respstd,'%4.3f') ];
    info.Temperature=[ num2str(info.physio.Temperature,'%4.2f') ' +- ' num2str(info.physio.Temperature,'%4.2f') ];
else
    info.resp=0; info.respstd=0; info.CO2=0; info.CO2std=0;
    info.CO2exp=0; info.CO2expstd=0;info.bpm=0;info.bpmstd=0;info.physio.Temperature=0;
end




function [out out2]=framestim(frame, stim)
%trouve le frame qui arrive immédiatement après la stimulation
%- On peut prendre un frame intermédiaire qui n'a pas été sauvegardé
%le dernier index dit après quel stim se trouve le frame
%out: même chose que frame mais avec un vecteur à al fin qui dit à quel stimulation on est rendu
%out2: les frame arrivant juste après la stimulation et son index initial
out=frame;
% for i1=1:
out2=zeros(size(stim,1),6);
out2=[];out3=[];
for i1=1:size(stim,1)
    a=find(frame(:,2)>stim(i1,2),1);
    if ~isempty(a);
        %     out2(i1,1:6)=[frame(a,1:5) a];
        dt=floor((frame(a,2)-stim(i1,2))*30);
        out2=[out2;frame(a,1:5)  a];
        out3=[out3;[frame(a,1:5) a] + [-dt -dt/30 -dt/30 0 0 0  ]];
        out(a:end,6)=i1;
    end
end
if isempty(out2)
    a=find(frame(:,2)>frame(1,2),1);
    out2=[frame(a,1:5) a];
end
if size(out2,1)<2
    b=find(frame(:,2)>(frame(end,2)-3),1);
    out2=[out2;frame(b,1:5) b];
end
