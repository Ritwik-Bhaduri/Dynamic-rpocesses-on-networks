G;
n = size(G);
n = n(1);
mu = 0;
lambda = 1;
state = datasample([0,1], n,'Weights',[0.99,0.01]);
% state_new = gillespie(state, G, mu, lambda);
% state_new = state_new(1:n)


iter = 4000;
repeat = 1;
t = zeros(iter,1);
prop_inf = zeros(iter,1);
    
hold on
for(j = 1:repeat)
    state_temp = state;
    for(i = 1:iter);
        state_temp = gillespie(state_temp(1:n), G, mu, lambda);
        prop_inf(i) = sum(state_temp(1:n))/n;
        t(i) = state_temp(n+1);
    end
    t = cumsum(t);
    Color=[1,0.5,0.5,0.3];
    plot(t,prop_inf,'color',Color);
end


function S = gillespie(state, G,  mu, lambda)
    n = size(G);
    n = n(1);
    S = state;
    propensity = zeros(n,1);    %propensity is the array of the transition rates of all the vertices
    for(i = 1:n)
        if(state(i)== 1)    propensity(i) = mu;
        else
            m = sum(G(i,:));
            propensity(i) = lambda * m;
        end
    end
    
    total_transition_rate = sum(propensity);
    rand1 = rand(1);
    t = -log(rand1)/total_transition_rate;      %t is the time of the first transition
    
    propensity = propensity/total_transition_rate;  %we rescale propensity array to get the vertex which changes
    temp = cumsum(propensity);
    rand2 = rand(1);
    flag = find(temp-rand2 > 0,1 ,'first');     %flag denotes index of the vertex that undergoes transition
    S(flag) = 1-S(flag);
    S = [S,t]; 
end