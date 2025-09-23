
clear


M = 50;
lambda = (1:(M-1)) / M;
nsample = 500;
n_cand = [1600 3200];
p = [1200];
b = 2.2;

rng(1)

tic

dnorm_cand = linspace(2,4,21);
dnorm_cand(1) = [];
delta_cand = 1-sqrt(dnorm_cand/4);
Delta = dnorm_cand([3 8 13]);

dnorm = zeros(length(delta_cand), 1);
size_zipf = zeros(length(Delta), length(delta_cand), length(n_cand));

for nn = 1:length(n_cand)
    n = n_cand(nn);
    N = floor((1 ./ rand(n,1)).^(1/(b-1)));

    lambda_thresholds = floor(lambda * floor(n/2));
    Nodd = N(1:2:(n-1));
    Neven = N(2:2:n);

    for dd = 1:length(delta_cand)
        delta = delta_cand(dd);
        mu = [(2-delta)*ones(p/2, 1); delta*ones(p/2, 1)]/p;
        mu0 = [delta*ones(p/2, 1); (2-delta)*ones(p/2, 1)]/p;


        dnorm(dd) = norm(mu-mu0)^2*p;
        [dd nn]




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

            rule(i,:) = (temp(end) > (Rn*10.2 + floor(n/2)*(floor(n/2)-1)/2*Delta));

        end
        size_zipf(:, dd, nn) = sum(rule) / nsample;

    end
end

toc


h = figure('Position', [100 100 500 300]);
h1 = plot(dnorm_cand,size_zipf(:,:,1),'k--');
hold on
h2 = plot(dnorm_cand,size_zipf(:,:,2),'k-');
hold off
line([Delta(1) Delta(1)], [0 0.2],'LineStyle', '--', 'Color', 'k')
line([Delta(2) Delta(2)], [0 0.2],'LineStyle', '--', 'Color', 'k')
line([Delta(3) Delta(3)], [0 0.2],'LineStyle', '--', 'Color', 'k')
set(h,'Units','Inches');
yline(0.05,'--')
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
legend([h1(1) h2(1)],{['n=',num2str(n_cand(1))],['n=',num2str(n_cand(2))]}, 'Location', 'northwest',FontSize=12);
print(h,['Setting1_b',num2str(10*b),'p',num2str(p)],'-dpdf','-r0')
box off




