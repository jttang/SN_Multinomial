
%% Zipf + Dirichlet

clear


rng(11)
M = 50;
lambda = (1:M) / M;
nsample = 50;
n = [1200];
p = [1200];
a_cand = 2.8;
delta_cand = 1;

size_zipf = zeros(length(delta_cand), length(a_cand));

tic


[i_idx, j_idx] = find(tril(ones(n), -2));
valid_pairs = j_idx <= floor((i_idx - 1) / 2);
i_idx = i_idx(valid_pairs);
j_idx = j_idx(valid_pairs);

for aa = 1:length(a_cand)
    a = a_cand(aa);
    N = floor((1 ./ rand(n,1)).^(1/(a-1)));
    Ni = N(i_idx);
    N2j_1 = N(2 * j_idx - 1);
    N2j = N(2 * j_idx);
    for dd = 1:length(delta_cand)

        dd

        delta = delta_cand(dd);

        mu = [(2-delta)*ones(p/3, 1); ones(p/3, 1); delta*ones(p/3, 1)]/p;

        mu0 = [delta*ones(p/3, 1); ones(p/3, 1); (2-delta)*ones(p/3, 1)]/p;

        Tn = zeros(nsample, 1);

        rule = zeros(nsample, 1);

        parfor i = 1:nsample
            % [dd aa i]
            X = mnrnd(N, mu');

            X_avg = (1 + sqrt((Ni .* N2j_1 + Ni .* N2j - 1) .* N2j ./ N2j_1)) .* X(2 * j_idx - 1, :) + ...
                (1 - sqrt((Ni .* N2j + Ni .* N2j_1 - 1) .* N2j_1 ./ N2j)) .* X(2 * j_idx, :);
            X_avg = X_avg ./ (N2j_1 + N2j);

            inner_term = sum((X(i_idx, :) ./ N(i_idx) - mu0') .* (X_avg - mu0'), 2);
            lambda_thresholds = floor(lambda * n);
            Sn = arrayfun(@(lambda) sum(inner_term(i_idx <= lambda)), lambda_thresholds)';

            floor_terms = floor(((3:n)-1) / 2);
            cumsum_floor_terms = cumsum(floor_terms);
            ratio = cumsum_floor_terms(lambda_thresholds-2);
            ratio = ratio' / ratio(end);

            Vn = sqrt(mean((Sn - ratio.* Sn(end)).^2));
            Tn(i) = Sn(end)/ Vn;

            Xtemp = sum(X);
            Ntemp = sum(N);
            expected = Ntemp .* mu0';
            hatmu = Xtemp / Ntemp;
            chi2_stat = sum((Xtemp - expected).^2 ./ expected);
            DRC = chi2_stat+1-sum(hatmu./mu0');
            rule(i) = (DRC > chi2inv(0.95, p-1)-(p-1));

        end
        size_zipf(dd, aa) = sum(Tn>5.88) / nsample;

    end
end


toc


size_zipf



