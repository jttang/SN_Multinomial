
%% Zipf + Dirichlet

clear

rng(11)
M_cand = [10 50 100];

nsample = 500;
n = [1600];
p = [1800];
a = 2.2;

tic

dnorm_cand = linspace(20/15,40/15,21);
dnorm_cand(1) = [];
delta_cand = 1-sqrt(dnorm_cand/8*3);
Delta = dnorm_cand([2 5 8]);
dnorm = zeros(length(delta_cand), 1);
size_zipf = zeros(length(Delta), length(delta_cand), length(M_cand));

for MM = 1:length(M_cand)
    M = M_cand(MM);

    lambda = (1:(M-1)) / M;
    N = floor((1 ./ rand(n,1)).^(1/(b-1)));

    lambda_thresholds = floor(lambda * floor(n/2));
    Nodd = N(1:2:(n-1));
    Neven = N(2:2:n);

    for dd = 1:length(delta_cand)
        delta = delta_cand(dd);
        mu = [(2-delta)*ones(p/3, 1); ones(p/3,1); delta*ones(p/3, 1)]/p;
        mu0 = [delta*ones(p/3, 1); ones(p/3,1); (2-delta)*ones(p/3, 1)]/p;


        dnorm(dd) = norm(mu-mu0)^2*p;
        [dd MM]

        rule = zeros(nsample, length(Delta));
        parfor i = 1:nsample
            X = mnrnd(N, mu');
            phi = (1 - sqrt((Nodd+Neven-1) .* Neven ./ Nodd)) ./ (Nodd+Neven);
            psi = (1 + sqrt((Nodd+Neven-1) .* Nodd ./ Neven)) ./ (Nodd+Neven);

            Y = phi .* X(1:2:(n-1),:) + psi .* X(2:2:n,:) - mu0';
            temp = p*sum(cumsum(Y).^2 - cumsum(Y.^2), 2) / 2;
            Tn = temp(lambda_thresholds);

            numerator = floor(lambda* floor(n/2)) .* (floor(lambda* floor(n/2))-1);
            ratio = numerator' / floor(n/2) / (floor(n/2-1));

            Rn = sqrt(sum((Tn - ratio.* temp(end)).^2)/(M-1));
            if M==10
                rule(i,:) = (temp(end) > (Rn*10.73 + floor(n/2)*(floor(n/2)-1)/2*Delta));
            elseif M==50
                rule(i,:) = (temp(end) > (Rn*10.20 + floor(n/2)*(floor(n/2)-1)/2*Delta));
            elseif M==100
                rule(i,:) = (temp(end) > (Rn*10.13 + floor(n/2)*(floor(n/2)-1)/2*Delta));
            end
        end
        size_zipf(:, dd, MM) = sum(rule) / nsample;
    end
end




size_zipf





