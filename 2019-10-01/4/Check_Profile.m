clear
load('1pipetteN100S100B2\winner300.mat')
fp_t = [];
n_t = [];
b_t = [];
for i = 1:100
    fp_t = [fp_t; fp{i}];
    n_t = [n_t; n{i}];
    b_t = [b_t; b{i}];
end
[fp_v, idx] = sort(fp_t, 'ascend');
n_v = n_t(idx);
b_v = b_t(idx);
bsize = 2e-4;
[fp_coarse1, edges] = histcounts(fp_v, (0:bsize:0.2)+bsize/2);
ctrs = edges(1:end-1)+bsize;
ctrs = ctrs(fp_coarse1>0);
fp_coarse = fp_coarse1(fp_coarse1>0);
fp_idx = cumsum([0 fp_coarse]);
fp_count = length(fp_idx)-1;
B_c = zeros(fp_count, 1);
fp_c = zeros(fp_count, 1);
for i = 1:fp_count
    B_c(i) = sum(b_v(fp_idx(i)+1:fp_idx(i+1)) .*n_v(fp_idx(i)+1:fp_idx(i+1)));
    fp_c(i) = ctrs(i);
end
%%
plot(fp_c, log(B_c), '.')
hold on
plot((0.02:0.01:0.2),-((0.02:0.01:0.2)-0.125).^2/1e-4+10)
hold off
ylim([0 15])