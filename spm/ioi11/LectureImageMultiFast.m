function [Imout frameout frameReadOut fileNo] = LectureImageMultiFast(Path,file,frameRead)
%Path :le repertoire contenant les fichiers .bin
% file: name of file   (ex: image)
%frameRead: first put FrameRead to -1, this return the first image and  frameReadOut)
% then take the proper line in FrameReadOut to select the good frame to read
% note, if frameRead contain 1 line, Imout est une matrix
%if frameRead is a vector, Imout is a cell
% frameReadOut : frames present in the file (frameNo  indexInFileNo fileNo)
% fileNo : number of the file containing image
persistent zzz filenameOld  filename00Old s1 s2 
if nargin==0
    clear zzz filenameOld  filename00Old s1 s2
    return
end
if ~strcmp(Path(end),'\')
    Path=[Path '\'];
end
if nargout>3
    fichier=dir([Path file '*.bin']);
    for i1=1:length(fichier)
        fileNo(i1)=str2num(fichier(i1).name(length(file)+1:end-4));
    end
    fileNo=sort(fileNo);
end
% oopen first file to see the size of the image
if ~strcmp([Path file],filename00Old) || frameRead(1)==-1
    fidA = fopen([Path file '0.bin']);
    I = fread(fidA,5,'int32');
    s1= I(3);
    s2= I(2);
    fclose(fidA);
end
if frameRead(1,1)==-1
    % extract frameReadOut i.e. all the frame present in the directory
    frameReadOut=zeros(100000,3);
    indA=1;
    for i1=1:length(fichier)
        fileNo(i1)=str2num(fichier(i1).name(length(file)+1:end-4));
        fidA = fopen([Path file   num2str(fileNo(i1)) '.bin']);
        EOF=0;
        indB=1;
        while  EOF==0
            try
                [ frameReadOut(indA,1)  count]= fread(fidA,1,'int32');
                [dummy]=fread(fidA,s1*s2+4,'int16');
                frameReadOut(indA,2)=indB;
                frameReadOut(indA,3)=fileNo(i1);
                indA=indA+1; indB=indB+1;
            catch
                EOF=1;
            end
        end
        fclose(fidA);
    end
    frameReadOut(indA:end,:)=[];
end

if frameRead(1,1)==-1
    frameRead = frameReadOut(1,:);
end

filename00Old=[Path file];
% lecture du .bin
for i1=1:size(frameRead,1)
    filename=[Path file num2str(frameRead(i1,3)) '.bin'];
    if ~strcmp(filename,filenameOld) %lecture of .bin if not already in memory
        fidA = fopen(filename);
        if fidA==-1
            disp(['fichier non présent:' filename])
            zzz=int16(zeros(10));
        else
            zzz=[];
            zzz=int16((fread(fidA,'int16')));
            zzz=reshape(zzz,(s1*s2+6),round(length(zzz)/(s1*s2+6)));            
            fclose(fidA);
        end
    else
        fidA=1; %si le fichier est déjà en mémoire on met le marqueur à 1 pour continuer l'analyse
    end
    
    if fidA~=-1   
        %disp(zzz(1,frameRead(i1,2)))
        frameout(i1)= frameRead(i1,1);
        Imout=double(reshape(zzz(7:end,frameRead(i1,2)),[s1 s2])'); %note :on tourne la figure
        
        if size(frameRead,1)>1
            Imoutcell{i1}= Imout;
        end
    else
        try 
            Disp2([],['frame : ' num2str(frameRead(i1)) ' non trouvé' ]) 
        catch
            disp([],['frame : ' num2str(frameRead(i1)) ' non trouvé' ]) 
        end
        Imout=zeros(s2, s1);
        frameout=0;
    end
    filenameOld=filename;
end

if size(frameRead,1)>1
    Imout= Imoutcell;
end

