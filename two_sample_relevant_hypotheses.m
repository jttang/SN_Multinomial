clear

rng(11)
M = 50;
lambda = (1:(M-1)) / M;
nsample = 500;
n_cand = [3200 6400];
p = [1200];
a = 2.2;
gamma = 1/2;

tic

dnorm_cand = linspace(20/15,40/15,21);
dnorm_cand(1) = [];
delta_cand = 1-sqrt(dnorm_cand/8*3);
Delta = dnorm_cand([2 5 8]);

dnorm = zeros(length(delta_cand), 1);
size_zipf = zeros(length(Delta), length(delta_cand), length(n_cand));

for nn = 1:length(n_cand)
    n = n_cand(nn);
    nX = n*gamma;
    nY = n - nX;
    N = floor((1 ./ rand(n,1)).^(1/(a-1)));
    Nodd = N(1:2:(n-1));
    Neven = N(2:2:n);
    NX = N(1:(n/2));
    NY = N((n/2+1):end);

    phi = (1 - sqrt((Nodd+Neven-1) .* Neven ./ Nodd)) ./ (Nodd+Neven);
    psi = (1 + sqrt((Nodd+Neven-1) .* Nodd ./ Neven)) ./ (Nodd+Neven);

    phiX = phi(1:(n/4));
    phiY = phi((n/4+1):end);
    psiX = psi(1:(n/4));
    psiY = psi((n/4+1):end);
    lambda_thresholds = floor(lambda * floor(nX/2));

    for dd = 1:length(delta_cand)
        delta = delta_cand(dd);
        muX = [(2-delta)*ones(p/3, 1); ones(p/3,1); delta*ones(p/3, 1)]/p;
        muY = [delta*ones(p/3, 1); ones(p/3,1); (2-delta)*ones(p/3, 1)]/p;


        dnorm(dd) = norm(muX-muY)^2*p;
        [dd nn]

        rule = zeros(nsample, length(Delta));
        parfor i = 1:nsample
            X = mnrnd(NX, muX');
            Y = mnrnd(NY, muY');

            Z = phiX .* X(1:2:(nX-1),:) + psiX .* X(2:2:nX,:) ...
                - phiY .* Y(1:2:(nY-1),:) - psiY .* Y(2:2:nY,:);
            temp = p*sum(cumsum(Z).^2 - cumsum(Z.^2), 2) / 2;
            Tn = temp(lambda_thresholds);

            numerator = floor(lambda* floor(n/4)) .* (floor(lambda* floor(n/4))-1);
            ratio = numerator' / floor(n/4) / (floor(n/4-1));

            Rn = sqrt(sum((Tn - ratio.* temp(end)).^2)/(M-1));

            rule(i,:) = (temp(end) > (Rn*10.20 + floor(n/4)*(floor(n/4)-1)/2*Delta));

        end
        size_zipf(:, dd, nn) = sum(rule) / nsample;
    end
end


h = figure('Position', [100 100 500 300]);
h1 = plot(dnorm_cand,size_zipf(:,:,1),'k--');
hold on
h2 = plot(dnorm_cand,size_zipf(:,:,2),'k-');
hold off
line([Delta(1) Delta(1)], [0 0.2],'LineStyle', '-.', 'Color', 'k')
line([Delta(2) Delta(2)], [0 0.2],'LineStyle', '-.', 'Color', 'k')
line([Delta(3) Delta(3)], [0 0.2],'LineStyle', '-.', 'Color', 'k')
yline(0.05,'-.')
set(h,'Units','Inches');
pos = get(h,'Position');
xlim([dnorm_cand(1) dnorm_cand(end)])
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
legend([h1(1) h2(1)],{'n=1600','n=3200'}, 'Location', 'northwest',FontSize=12);
print(h,['1two_b',num2str(10*a),'p',num2str(p)],'-dpdf','-r0')
box off

