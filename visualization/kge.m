function k = kge(obs, sim)
% Define Kling-Gupta Efficiency (KGE) function
r = corr(obs, sim, 'rows', 'complete');
alpha = std(sim) / std(obs);
beta = mean(sim) / mean(obs);
k = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
end
