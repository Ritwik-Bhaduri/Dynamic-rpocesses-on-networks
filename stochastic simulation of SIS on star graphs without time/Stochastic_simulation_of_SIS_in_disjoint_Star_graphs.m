% state represents the array of states of all the vertices in entire graph
% state is represented in 0 and 1; 0 for susceptible and 1 for infected
% There are B many stars in the graph In avertex the first vertex is the
% cetral node and the leaves come next

%% Initialisation:
B = 50;
state = zeros(B * 5,1);

mu = 1;
lambda = 5;

state = datasample([0,1], (B * 5));
sum(state) / (B * 5)

S = zeros(4,5); % vertices in S(i,j) have i neighbours and (m-1) neighbours
I = zeros(4,5); % vertices in I(i,j) have i neighbours and (m-1) neighbours

adjacency_matrix = zeros(5,5);
adjacency_matrix(1,:) = [0,1,1,1,1];
adjacency_matrix(:,1) = [0,1,1,1,1];
A = adjacency_matrix;
A = repmat(A,B,1);

%% Counting the ratios (counts)
% for(i = 1:B)
%     local_state = state(5 * (i-1) + 1 : 5 * i ).';
%     local_state = reshape(local_state, 1, 5);
%     for(j = 1:5)
%         v = A(5 * (i-1) + j, :);
%         n_inf_ngbr = sum(local_state .* v);
%         n_sus_ngbr = sum(v) - n_inf_ngbr;
%         k = sum(v);
%         m = n_inf_ngbr;
%         if(local_state(j) == 1) I(k,m + 1) = I(k,m + 1) + 1;
%         else S(k,m + 1) = S(k,m + 1) + 1;
%         end
%     end
% end
% s = S / (B * 5);
% i = I / (B * 5);
% S;
% I;

%% Simulation

mu
lambda
iter = 1000000;
prop_inf = zeros(floor(iter),1);

for(i = 1:iter)
    state = simulate2(state, A, mu, lambda);
    prop_inf(i) = sum(state) / (B * 5);
    
    i = i + 1;
end
plot(prop_inf)

% proportion_of_infected_vertices = sum(state) / (B * 5)
% proportion_of_succeptible_vertices = 1 - proportion_of_infected_vertices
