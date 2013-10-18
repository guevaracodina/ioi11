function [L,H,F,C,T,B,R,tf,codage,file] = ioi_read_blk_info(file, affiche)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonction [L,H,F,C,T,B,R,tf,codage,file] = ioi_read_blk_info(file, affiche);
% Lecture des informations dans les fichiers sources de type .BLK et
% affichage optionnel de ces paramètres
%
% Paramètres d'entrée :
% file = fichier source de type .BLK
% affiche = Affichage (1) ou non (0) des différents paramètres de sortie
%           dans le Command Window de Matlab
%
% Paramètres de sortie :
% L = largeur des images
% H = hauteur des images
% F = nombre de frame par condition
% C = nombre de condition
% T = nombre de trial par block
% B = nombre de block (dans le repertoire)
% R = nombre de repetition
% tf = duree d'un segment (stim + interstim)
% codage = codage du fichier ('uint16' = 16 bits (2 octets), 'uint32' = 32 bits (4
% octets)
%
% Affichage :
% Affichage des paramètres de sortie si désiré dans le Command Window de
% Matlab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Si le nombre d'argument est de 1, on affiche pas les paramètres de sortie
if nargin==1
    affiche=0;
end

%Lecture de différents paramètres de sorties
fid = fopen(file);
A = fread(fid,500,'uint16');
L = A(19);
H = A(21);
F = A(23);
C = A(25);
T = A(467);
tf = A(165);
codage = A(17);
fclose(fid);
courant = pwd;

if sum(file == '\') ~= 0
    travail = file(1:max(find(file == '\')));
    vrai_file = file(max(find(file == '\'))+1:length(file));
else
    travail = pwd;
    vrai_file = file;
end;

%Détermination du nombre de bloc dans le répertoire
cd(travail)
rep = dir;
[nf,~] = size(rep);
B = 0;
prefix1 = vrai_file(1:length(vrai_file)-5);
for i = 1:nf
    [~,nn] = size(rep(i).name);
    if nn == length(vrai_file) || nn == length(vrai_file) + 1 
        prefix2 = rep(i).name(1:length(vrai_file)-5);
        if prefix1 == prefix2
            B = B + 1;
        end;
    end;
end;
cd(courant)
R = B * T;

%Affichage optionnel des paramètres de sortie
if affiche==1;
disp(['noms des fichiers = ' vrai_file])
disp(['largeur = ', num2str(L)])
disp(['hauteur = ', num2str(H)])
disp(['nombre de frames = ', num2str(F)])
disp(['nombre de conditions = ', num2str(C)])
disp(['nombre de trials par block = ', num2str(T)])
disp(['nombre de blocks = ', num2str(B)])
disp(['nombre de repetitions = ', num2str(R)])
disp(['codage (octets) = ', num2str(codage)])
end;
