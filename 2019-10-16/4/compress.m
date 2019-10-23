function [fp, b, n, n_genos] = compress(fp, b, n, dFP_dL)

% cw_params = [compressWinner, digits_fp, digits_b];
% dFP_dL = [digits_fp, digits_b]
% 'compresses' the three arrrays by binning fp and b.
nbins_fp = 10.^dFP_dL(1);
nbins_b = 10.^dFP_dL(2);
lossy_fp = round(nbins_fp * fp) / nbins_fp;
lossy_b = round(nbins_b * b) / nbins_b;
FLN = sortrows([lossy_fp,lossy_b,n],[1,2]);
for i = 2 : size(FLN, 1)%length(FLN)
    if isequal(FLN(i,1:2),FLN(i-1,1:2))
        FLN(i,3) = FLN(i,3) + FLN(i-1,3);
        FLN(i-1,3) = 0;
    end
end
FLN(FLN(:,3) == 0,:) = [];
fp = FLN(:,1);
b = FLN(:,2);
n = FLN(:,3);
n_genos = length(fp);
end