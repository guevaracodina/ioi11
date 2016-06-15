%% Read Anatomical Data
dataName = 'C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Data';
resultsName = 'C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Data\Atlas';
subjectName{1} = 'FB25E02';
subjectName{2} = 'FB25E03';
subjectName{3} = 'FB26E03';
subjectName{4} = 'FB26E04';
subjectName{5} = 'FB26E05';
subjectName{6} = 'FB26E06';

%% Anatomical
% Bregma est environ à la ligne 250, colonne 400 et Lambda: 250, 150.
% pseudo_anat = mat2gray(squeeze(mean(HbO,1)));
for iSubjects=1:numel(subjectName)
    h = open(fullfile(dataName,[subjectName{iSubjects} '.fig']));
    set(h,'units','inch')
    anat=getimage;
    imwrite(mat2gray(anat), fullfile(resultsName,[subjectName{iSubjects} '_anat_S01.png']),...
        'BitDepth',16,'XResolution',300,'YResolution',300)
    % ,'XResolution',300,'YResolution',300
    close(h)
end
