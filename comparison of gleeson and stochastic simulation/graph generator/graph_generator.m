%% Erdos renyi graph generator

n = 5000;
p = 0.02;
G = rand(n,n) < p;
G = triu(G,1);
G = G + G';
plot(graph(G));

%% Configuration model generator
% the choices for degree distributions are : 'uniform','normal','binomial','exponential'

n = 100;
dist = 'exponential'
G = randomGraphDegreeDist(n, dist);
plot(graph(G));

%% Complete graph

n = 70;
G = zeros(n , n);
for(i = 1:n)
    for(j = 1:n)
        if(i ~= j) G(i,j) = 1;
        end
    end
end
plot(graph(G));

%% Configuration model for powerlaw distribtion
%the probabilities are proportional to p^(-alpha)

n = 1000;
alpha = 2.5;
c = 0;
lower_limit = 3;
upper_limit = 20;
values = [lower_limit : upper_limit];
probability = zeros(upper_limit-lower_limit+1,1);
for(i = values)
    i;
    probability(i - lower_limit + 1) = i^(-alpha);
end
probability = probability/sum(probability);


Nseq = zeros(length(values));
while not(isGraphic(Nseq)); Nseq = weightedRandomSample(n,values,probability); end
G = randomGraphFromDegreeSequence(Nseq);
plot(graph(G))