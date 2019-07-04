%simulate2 simulates SIS dynamics for a collection of star graphs
%   simulate2(state, A, mu, lambda) simulates SIS dynamics on a collection
%   of star graphs stochastically using a variant of Gillespie's algorithm. 
%   A star graph is a central node with 4 leaf nodes.
%   The function will simulate only one step of the process.

function S = simulate2(state, A, mu, lambda)
    B = length(state)/5;
    S = state;
    c = datasample([1:B],1); 
    local_state = state((5 * (c-1) + 1) : (5 * c)).';
    local_state = reshape(local_state,1,5);
    local_A = A(((5 * (c-1) + 1) : (5 * c)),:);
    S((5 * (c-1) + 1) : (5 * c)) = simulate1(local_state, local_A, mu, lambda);
end