function imageOut=ioi_read_blk_images(file,bz,bxy,r,base,filtre,FM,N,cible,sens,pourcent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonction []=Compile_Retinotopie_Dominance_oculaire(file,bz,bxy,r,base,filtre,FM,N,cible,sens,pourcent)
% Compile les fichier de stimulation continue de VDAQ (.BLK) en .MAT et
% effectue des opérations de base comme le binning temporel, découpage
% spatial et temporel, correction d'illumination et filtrage spatial.
%
% Paramètres d'entrée :
% file = fichier source (BLK)
% bz = decoupage temporel [MIN MAX], si [0] = pas de decoupage
% bxy = decoupage spatial [xMIN xMAX yMIN yMAX], si [0] = pas de decoupage
% r = binning temporel, si [1] = pas de binning
% base = soustraction de la ligne de base 
%            [N] = moyenne sur ''N'' pixel
%            [0] = pas de ligne de base
%            [-1] = soustraction de la moyenne
%            [-2] = soustraction instabilite de lampe sur toute l'image
%            [-3] = soustraction instabilite de lampe a partir d'une ROI 
%                   reference a cliquer
%                       - si possible, ne doit pas contenir d'activation
%                       - si possible, doit avoir une illimation similaire 
%                         a la zone activee
%                       - si possible, ne doit pas contenir de gros VS
%            [a b c d] = idem mais roi def
% filtre = filtrage spatial [coupure_moyen, spatial_gaussien, 
%          coupure_gaussien], si [0] = pas de filtrage
% FM = limite des fréquences conservées lors de la Transformée de Fourier
% N = nombre d'étapes de traitement
% cible = fichier .MAT cible
% sens = sens du signal (-1 = neg)
% pourcent = 1 si mise en pourcentage, 0 autrement
%
% Paramètres de sortie : aucun
%
% Affichage :
% Affichage de la progression de la conversion dans le Command Window de
% Matlab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ouverture
if strcmp(file(end-2:end),'BLK') == 1
    disp('- ouverture de fichier VDAQ')
    [L,H,F,~,~,imageOut,~,~,codage] = ioi_read_blk_info(file,0);
    fid = fopen(file,'r');
    fseek(fid,1716,-1);
    if codage == 2
        disp('>>>>>>codage en 16 bits')
        imageOut = single(fread(fid,(L*H*F),'uint16'));
    elseif codage == 4
        disp('>>>>>>codage en 32 bits')
        imageOut = single(fread(fid,(L*H*F),'uint32'));
    end;
    fclose(fid);
    imageOut=reshape(imageOut,L,H,F,1);
else
    disp('- ouverture de matrice matlab dans le workspace')
    load(file)
    [L,H,F] = size(imageOut);
end;

% Rotation des images en position anatomique
B2=zeros(H,L,F);
for i=1:F
    B2(:,:,i)=flipud(rot90(imageOut(:,:,i)));
end
imageOut=B2;
clear B2
[L,H,F]=size(imageOut);
%% decoupage et binning temporel
if sum(bz) ~= 0
    disp('- decoupage temporel')
    imageOut = imageOut(:,:,bz(1):bz(2));
    [L,H,F] = size(imageOut);
end;
if sum(bxy) ~= 0
    disp('- decoupage spatial')
    imageOut = imageOut(bxy(1):bxy(2),bxy(3):bxy(4),:);
    [L,H,F] = size(imageOut);
end;
if r ~= 1
    disp('- rebinning temporel')
    imageOut=reshape(imageOut,L*H,F,1);
    imageOut = imresize(imageOut,[L*H round(1/r*F)]);
    F = round(1/r*F);
    imageOut=reshape(imageOut,L,H,F,1);
end;
%% soustraction ligne de base (correction d'illumination)
if length(base) == 1
    if base > 0
        disp('- soustraction ligne de base classique')
        segment = L*H/N;
        tic
        imageOut=reshape(imageOut,L*H,F,1);
        for i = 1:N
            texte = ['etape ' num2str(i) ' sur ' num2str(N)];
            disp(texte)
            REF = (conv2(imageOut((i-1)*segment+1:segment*i,:),ones(1,base)/base,'same') ./ ...
                 conv2(ones(size(imageOut((i-1)*segment+1:segment*i,:))),ones(1,base)/base,'same'));
             if pourcent == 1
                 disp('     /// mise en pourcentage')
                 imageOut((i-1)*segment+1:segment*i,:) = (imageOut((i-1)*segment+1:segment*i,:) - REF) ./REF;
             else
                 disp('     ///soustraction ligne de base')
                 imageOut((i-1)*segment+1:segment*i,:) = imageOut((i-1)*segment+1:segment*i,:) - REF;
             end;
        end;
        imageOut=reshape(imageOut,L,H,F);
        toc
    elseif base == -1 
        disp('- soustraction ligne de base de la moyenne')
        m = mean(imageOut,3);
        for i = 1:F
            imageOut(:,:,i) = (imageOut(:,:,i) - m)./m;
        end;
    elseif base == -2 
        disp('- soustraction instabilite lumiere a partir de toute l''image')
        M = mean(imageOut,3);
        C = zeros(size(imageOut));
        TT = squeeze(mean(mean(imageOut,1),2));
        TT = TT / mean(TT);
        for i = 1:L
            for j = 1:H
                C(i,j,:) = M(i,j) * reshape(TT,1,1,F); 
            end;
        end;
        imageOut = imageOut - C;
        clear C
    elseif base == -3
        disp('- soustraction instabilite lumiere, clicage')
        disp(' >>>>>  CLIQUER 2 POINTS POUR DELIMITER UNE ROI REFERENCE')
        imagesc(imageOut(:,:,1)), colormap jet, title('CLIQUER 2 POINTS POUR DELIMITER UNE ROI REFERENCE')        
        [x1,y1] = ginput(1); x1 = round(x1); y1 = round(y1);
        [x2,y2] = ginput(1); x2 = round(x2); y2 = round(y2);
        close
        M = mean(imageOut,3);
        C = zeros(size(imageOut));
        TT = squeeze(mean(mean(imageOut(y1:y2,x1:x2,:),1),2));
        TT = TT / mean(TT);
        for i = 1:L
            for j = 1:H
                C(i,j,:) = M(i,j) * reshape(TT,1,1,F); 
            end;
        end;
        imageOut = imageOut - C;
        clear C
    end;
else
    disp('- soustraction instabilite lumiere dans une ROI')
    M = mean(imageOut,3);
    C = zeros(size(imageOut));
    TT = squeeze(mean(mean(imageOut(base(1):base(2),base(3):base(4),:),1),2));
    TT = TT / mean(TT);
    for i = 1:L
        for j = 1:H
            C(i,j,:) = M(i,j) * reshape(TT,1,1,F); 
        end;
    end;
    imageOut = imageOut - C;
    clear C
end;


%% filtrage 2D spatial classique par bout
if sum(filtre) ~= 0
    disp('- filtrage spatial')
    tic
    for i = 1:N
        imageOut(:,:,(i-1)*(F/N)+1:i*(F/N)) = OIA_pttt('data',...
            imageOut(:,:,(i-1)*(F/N)+1:i*(F/N)),F/N, 'filtrage',filtre);
    end;
    toc
end;
%% FFT sequencielle dans la dimension des frames
% disp('- transformee de fourier')
% imageOut=reshape(imageOut,L*H,F,1);
% [X,~] = size(imageOut);
% segment = X/N; 
% A = zeros(X,FM,'single');
% phi = zeros(X,FM,'single');
% disp(['--------nombre de segment : ',num2str(N)])
% for i = 1:N
%     tic
%     disp(['segment ',num2str(i)])
%     TF = zeros(round(segment),F);
%     TF = fft(sens*imageOut(round((i-1)*segment+1):round(i*segment),:),F,2);
%     A(round((i-1)*segment+1):round(i*segment),:) = abs(TF(:,1:FM));
%     phi(round((i-1)*segment+1):round(i*segment),:) = angle(TF(:,1:FM));
%     disp(['     -> (s) ',num2str(toc)])
% end;
% A = reshape(A,L,H,FM) / F;
% phi = reshape(phi,L,H,FM);
% imageOut=reshape(imageOut,L,H,F);
%% Enregistrement
% save(cible,'A','phi','L','H','F','file')
