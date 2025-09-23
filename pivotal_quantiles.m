
%% W_1
clear
clc

rng(1)

M = 10;     % = 20; 50; 100;
n = 10000;
lambda = (1:M) / M; 
Btemp = randn(n, M);
temp = cumsum(Btemp.*sqrt(diff([0 lambda].^2)), 2);
temp2 = (temp - temp(:, end) * lambda.^2).^2;

W1 = temp(:, end) ./ sqrt(mean(temp2, 2));

round(quantile(W1,[0.99, 0.975, 0.95, 0.90]),5)


%%

rng(1)

lambda = (1:M) / M; 
Btemp = randn(n, M);
temp = cumsum(Btemp.*sqrt(diff([0 lambda])), 2);
temp2 = (temp - temp(:, end) * lambda).^2 .* lambda.^2;
W2 = temp(:, end) ./ sqrt(mean(temp2, 2));

round(quantile(W2,[0.99, 0.975, 0.95, 0.90]),5)