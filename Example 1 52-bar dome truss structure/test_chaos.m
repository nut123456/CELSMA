% Generate 500 values using the Logistic map with a scaling factor of 10
Max_Iteration = 500;
Value = 10;
G = zeros(1, Max_Iteration); % Preallocate G for efficiency

for curr_iter = 1:Max_Iteration
    G(curr_iter) = chaos(5, curr_iter, Max_Iteration, Value);
end

% Plot the results
plot(G);
title('Logistic Map Chaotic Sequence');
xlabel('Iteration');
ylabel('Scaled Value');
