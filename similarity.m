function [similar] = similarity(theta)

N = length(theta(:,1));

% Get magnitude of estimate theta for each model
for i = 1:N
    magM(i) = sqrt(sum(theta(i,:).^2));
end

% Dot product
for i = 1:N
    for j = 1:N
        dot_prod(i,j) = theta(i,:)*theta(j,:)';
        Mg(i,j) = magM(i)*magM(j);
    end
end

% Model Cosine Similarity
similar = triu(dot_prod./Mg);


