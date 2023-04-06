function [D] = inverse_D22(Pk, Hn, n)


% Get H2
Hn1 = Hn(:,1);

% Get invH2*Pk*H2 
D = inv(Hn1'*Pk*Hn1);

% Recursive update
for j = 2:n

    % Hn-1 and hn
    Hn1 = Hn(:, 1:j-1);
    hn = Hn(:, j);

    % Common terms
    hnPk = hn'*Pk; 
    % ( 1 x T )

    Hn1Pkhn = Hn1'*hnPk';
    % ( n-1 x  1 )

    temp = D*Hn1Pkhn;
    % ( n-1  x  1 )

    % Get last element scalar c
    c = hnPk*hn - Hn1Pkhn'*temp;
    c = 1/c;

    % Get b vector
    b = -temp*c;

    % Get A
    A = D + c*(temp*temp');

    % Update new D
    D = [A, b; b', c];

end



end