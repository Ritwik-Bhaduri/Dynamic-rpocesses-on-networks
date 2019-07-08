%% Initialisation:
B = 1; % B is the number of star graphs
iter = 100000; % iter gives the number of processes that are simulated for a particular initial state.
mu = 1;
lambda = 4.5;

state = datasample([0,1], (B * 5));
initial_proportion_of_infected_vertices = sum(state) / (B * 5)

S = zeros(4,5); % vertices in S(i,j) have i neighbours and (m-1) neighbours
I = zeros(4,5); % vertices in I(i,j) have i neighbours and (m-1) neighbours

adjacency_matrix = zeros(5,5);
adjacency_matrix(1,:) = [0,1,1,1,1];
adjacency_matrix(:,1) = [0,1,1,1,1];
A = adjacency_matrix;
A = repmat(A,B,1);


%% Simulation
state = datasample([0,1], (B * 5));
t_max = 2 * 10 ^ 5;
t_com = 1:t_max;
prop_inf_com = zeros(t_max,1);
iterations = 5; 
% iterations gives the number of repeatations for averaging
for(k = 1 : iterations)
    state1 = state;
    prop_inf = zeros(iter,1);
    t = zeros(iter,1);
    for(i = 1:iter)
        arr = simulate2withtime(state1, A, mu, lambda);
        state1 = arr(1: B * 5);
        t(i) = arr(B * 5 + 1);
        prop_inf(i) = sum(state1) / (B * 5);
    end
    t_cum = cumsum(t);
    for(i = 1 : t_max)
        j = find(t_cum < i, 1, 'last');
        prop_inf_com(i) = prop_inf_com(i) + sum(prop_inf(j));
    end
end
prop_inf_com = prop_inf_com / iterations;
plot(t_com,prop_inf_com);
