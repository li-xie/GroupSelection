% use the symmetric mutation spectrum

function [fp_mut] = mut_spec(params,fp_background)
sp = params(1); % mean increase of a mutation
sn = params(2); % mean decrease of a mutation
fp_min = params(3); % epistatic factor
fp_max = params(4); % maximal value of fp
% sp = sp0 ./ (1 + g * (fp_background - 0.13) / 0.13);
% sn = sn0 * (1 + g * (fp_background - 0.13) / 0.13);
% the normalization constant of the mutation spectrum
nc = sp + sn;
% the division point of the spectrum
divalue = sn/nc;

u = rand(size(fp_background));
delta_fp = zeros(size(fp_background));

% inverse transformation
idx = find(u <= divalue);
if ~isempty(idx)
    delta_fp(idx) = sn *log(u(idx) *nc/sn);
end
idx = find(u > divalue);
if ~isempty(idx)
    delta_fp(idx)=-sp*log((1-u(idx)) *nc/sp);
end

fp_mut = fp_background + delta_fp;
fp_mut(fp_mut > fp_max) = fp_max;
fp_mut(fp_mut < fp_min) = fp_min;
end