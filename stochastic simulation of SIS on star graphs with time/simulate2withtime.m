function S = simulate2withtime(state, A, mu, lambda)
    B = length(state)/5;
    S = state;
    c = datasample([1:B],1); 
    local_state = state((5 * (c-1) + 1) : (5 * c)).';
    local_state = reshape(local_state,1,5);
    local_A = A(((5 * (c-1) + 1) : (5 * c)),:);
    arr = simulate1temp(local_state, local_A, mu, lambda);
    S((5 * (c-1) + 1) : (5 * c)) = arr(1:5);
    S = [S, arr(6)];
end
