%adjacency matrix for star graph taking the centre vertex as the first node
adjacency_matrix = zeros(5,5);
adjacency_matrix(1,:) = [0,1,1,1,1];
adjacency_matrix(:,1) = [0,1,1,1,1];
A = adjacency_matrix;

S0 = [1,0,0,0,0];    %initialization:  infected are 1 and succeptibles are 2
mu = 0.05;
lambda = 1;

%%
S = simulate1(S0, A, mu, lambda);
disp(S)
S = S(1:5)
for(i = 1:10000)
    S1 = simulate1(S, A, mu, lambda);
    disp(S1);
    S = S1(1:5);
    if(S == [0,0,0,0,0]) 
        disp(['step at which process stops is  ',num2str(i+1)]);
        break
    end
end

%%
function S = simulate1(state, A, mu, lambda)
    if(state == [0,0,0,0,0]) 
        disp('All the nodes are succeptible. So, no further change will happen')   
    else
        count = zeros(5,1);
        %count is the number of 2s in each row after multiplying by the State vector. 
        %This gives us the number of I-S edges for each vertex
        for(i = 1:5)
            count(i) = sum((repelem(1,5)-state) .* A(i,:));
        end
        Total_IS = state * count; %this gives the number of lambdas i.e. the number of Infection processes possible
        Total_mu = sum(state); %this gives the number of mu's i.e. the number of recoveries possible
        r1 = rand(1);
        if(r1 <= (Total_mu * mu)/(Total_mu * mu + Total_IS * lambda)) flag = 0;
        else    flag = 1; end      %flag is 0 if a recovery has occurred else if an infection has occurred the flag is 1
        if(flag == 0)
            v = 1:5;
            infected_vertices = v(state == 1);
            recovering_vertex = datasample(infected_vertices,1);
            state(recovering_vertex) = 0;
        else
            n_IS = state .* count.';
            partition = cumsum(n_IS)/max(cumsum(n_IS));
            r3 = rand(1);
            infecting_vertex = find((partition - r3) > 0, 1);
            v = (1:5) .* A(infecting_vertex,:) .* (repelem(1,5) - state);
            succeptible_vertices = v(v > 0);
            infected_vertex = datasample(succeptible_vertices,1);
            state(infected_vertex) = 1;
        end
    end
    S = state;
end