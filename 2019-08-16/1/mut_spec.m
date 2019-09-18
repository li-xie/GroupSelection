% use the symmetric mutation spectrum

function [fp_mut] = mut_spec(params,fp_background)
sp0 = params(1); % mean increase of a mutation
sn0 = params(2); % mean decrease of a mutation
g = params(3); % epistatic factor
fp_max = params(4); % maximal value of fp
sp = sp0 ./ (1 + g * (fp_background - 0.13) / 0.13);
sn = sn0 * (1 + g * (fp_background - 0.13) / 0.13);
% the normalization constant of the mutation spectrum
nc = sp + sn .* (1 - exp(-1 ./ sn));
% the division point of the spectrum
divalue = sn .* (1 - exp(-1 ./ sn)) ./ nc;

u = rand(size(fp_background));
delta_fp = zeros(size(fp_background));

% inverse transformation
idx = find(u <= divalue);
if ~isempty(idx)
    delta_fp(idx) = sn(idx) .* log(u(idx) .* (sp(idx) + sn(idx) .* ...
        (1 - exp(-1 ./ sn(idx)))) ./ sn(idx) + exp(-1 ./ sn(idx)));
end
idx = find(u > divalue);
if ~isempty(idx)
    delta_fp(idx)=-sp(idx).*log((sp(idx)+sn(idx).*(1-exp(-1./sn(idx))))./sp(idx).*(1-u(idx)));
end

fp_mut = fp_background .* (1 + delta_fp);
fp_mut(fp_mut > fp_max) = fp_max;
end