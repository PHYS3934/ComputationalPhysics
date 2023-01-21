function M = construct_update_matrix(L,whatMethod,g)
% Constructs an update matrix
%-------------------------------------------------------------------------------

switch whatMethod
case 'ftcs_advection'
    % Diagonal matrix to implement the centred spatial first-derivative
    D = diag(ones(L-1,1),+1) - diag(ones(L-1,1),-1);

    % Additional elements for periodic boundary conditions
    D(1,L) = -1;
    D(L,1) = +1;

    % Update matrix
    M = eye(L) - 0.5*g*D;

case 'lax_advection'
    % Diagonal matrix to implement the centred spatial first-derivative
    D = diag(ones(L-1,1),+1) - diag(ones(L-1,1),-1);

    % Additional elements for periodic boundary conditions
    D(1,L) = -1;
    D(L,1) = +1;

    % Lax averaging matrix
    A = 0.5*abs(D);

    % Update matrix
    M = A - 0.5*g*D;
otherwise
    error('Unknown method: ''%s''',whatMethod)
end

end
