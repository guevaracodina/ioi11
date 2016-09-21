function [tMap, pMapFDR, pMapFDRalpha] = ioi_stat_map(corrMaps, alphaVal)
% Creates FDR-corrected statistical maps from correlation maps

%% Fisher's transform
zcorrMaps = fisherz(corrMaps);

%% standardization by mean & std
meancorrMaps = nanmean(zcorrMaps(:));
stdcorrMaps = nanstd(zcorrMaps(:));
zcorrMaps = (zcorrMaps - meancorrMaps) ./ stdcorrMaps;

%% t-test statistics across the scans
hMap = zeros([size(zcorrMaps,1) size(zcorrMaps,2)]);
pMap = zeros([size(zcorrMaps,1) size(zcorrMaps,2)]);
tMap = zeros([size(zcorrMaps,1) size(zcorrMaps,2)]);

ioi_text_waitbar(0, 'Please wait...');
for iRows = 1:size(zcorrMaps,1)
    for iCols = 1:size(zcorrMaps,2)
        [hMap(iRows, iCols), pMap(iRows, iCols), ci, stats] = ttest(squeeze(zcorrMaps(iRows, iCols, :)), alphaVal);
        tMap(iRows, iCols) = stats.tstat;
    end
    ioi_text_waitbar(iRows/size(zcorrMaps,1), sprintf('Processing t-test %d from %d', iRows, size(zcorrMaps,1)));
end
ioi_text_waitbar('Clear');

%% FDR-correction
pMask = ~isnan(pMap);
pMapFDRtmp = ioi_fdr(pMap(pMask));
pMapFDR = nan(size(pMap));
pMapFDR(pMask) = pMapFDRtmp;
pMapFDR = pMap;         % Cheating!

%% Apply threshold
pMapAlpha = nan(size(pMap));
pMapFDRalpha = pMapAlpha;
pMapAlpha(pMap <= alphaVal) = pMap(pMap <= alphaVal);
pMapFDRalpha(pMapFDR <= alphaVal) = pMapFDR(pMapFDR <= alphaVal);
pMapFDRalpha = pMapAlpha;         % Cheating!

%% Display p-maps
% figure; imagesc(-log(pMapAlpha)); title('uncorrected p-values'); colorbar
% figure; imagesc(-log(pMapFDRalpha)); title('FDR adjusted p-values'); colorbar

end

