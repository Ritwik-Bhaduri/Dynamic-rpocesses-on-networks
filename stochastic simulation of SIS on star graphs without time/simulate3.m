%simulate3 Simulates accelerated SIS dynamics for a collection of star graphs
%   simulate3(state, A, mu, lambda) simulates SIS dynamics on a collection
%   of star graphs stochastically using a variant of Gillespie's algorithm.
%   A star graph is a central node with 4 leaf nodes.
%   In every timestep simulate3 simulates one process in each of the star
%   graphs as they are independent. In this way it arrives at the final
%   limitng state of the process faster than simulate2 
%   (which is the normal SIS dynamics).
%   (Note: The proportions of the intermeiate steps might not agree with that of pure SIS)
%   The function will simulate only one step of the process.

function S = simulate3(state, A, mu, lambda)
    B = length(state)/5;
    S = zeros(B * 5,1);
    for(i = 1:B)
        local_state = state((5 * (i-1) + 1) : (5 * i)).';
        local_state = reshape(local_state,1,5);
        local_A = A(((5 * (i-1) + 1) : (5 * i)),:);
        S((5 * (i-1) + 1) : (5 * i)) = simulate1(local_state, local_A, mu, lambda);
    end
    S = S;
end