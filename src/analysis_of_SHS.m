%% Plot a 3D magnitude graph to have a idea of the S'*S matrix

SHS = S_normalized'*S_normalized;

% Calculate the magnitudes of the complex entries
magnitudes = abs(SHS);

% Create a meshgrid for the matrix indices
[x, y] = meshgrid(1:size(SHS, 2), 1:size(SHS, 1));

% Create the 3D surface plot
figure;
surf(x, y, real(SHS), 'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('Column Index');
ylabel('Row Index');
zlabel('Magnitude');
title('3D Plot of Magnitudes in Complex Matrix');

% Customize the plot appearance (optional)
colormap('jet'); % You can choose a different colormap if you prefer
colorbar; % Add a colorbar to indicate magnitude values

% Rotate the plot to a better view (optional)
view(45, 30); % You can adjust the azimuth and elevation angles as needed

%% Result evaluation (x axis is sparsity)
k = 0:1:N-1;

nnze = zeros(N,1);
for i = 1:1:N
    % MSE(i) = (W(:,i) - w)'*S'*S*(W(:,i) - w)/L;
    nnze(i) = nnz(W(:,i));
end
subplot(4,1,1)
p1 = semilogy(k,flip(abs(MSE)));
% plot(k,flip(abs(MSE)))
xlabel('Sparsity')
ylabel('MSE')
title('The MSE of the threshold pattern vs sparsity')

hold on

p2 = semilogy(k,flip(abs(MSE_app)),"-.");

p3 = semilogy(k,flip(abs(MSE_new)),"--");

legend([p1,p2,p3], ["MSE", "Approximation", "New"]);

subplot(4,1,2)
stem(flip(MSE - MSE_new))
ylabel('The difference in magnitude ')
title('The difference between MSE calculated with (Sw - Sw_t) and S(w - w_t)')

subplot(4,1,3)
stem(flip(MSE_new - MSE_app))
ylabel('The difference in magnitude ')
title('The difference between MSE and its approximation')

subplot(4,1,4)
plot(k,flip(nnze))
xlabel('Sparsity')
ylabel('Number of non-zeros entries')
title('The number of nonzero element vs sparsity')

figure
plot(k,flip(100*abs(MSE - MSE_app)/MSE));
ylabel('Relative error in percentage')
title('The relative error of approximation of MSE')

%% Check the percentage of maginitude of diagonal entries occupy in the whole matrix
% sum of diagonal entries
sum_diag = trace(magnitudes);

% sum of all entries
sum_all = sum(magnitudes,"all");

percentage_diag_all = 100*sum_diag/sum_all

% This test is not a very good test

%% Check the goodness of representing S^H*S with only diagonal entries 
wtest = ones(N,1);
result_all = wtest'*SHS*wtest;

result_diag = wtest'*diag(diag(SHS))*wtest;

percentage_diag_all = 100*result_diag/result_all


%% Check the difference between desired pattern and the generated pattern
figure
stem(desired_pattern_modified - abs(S_newn*wn))
title('The difference between desired pattern and generated pattern at L dimensions')
% The most of the epsilon (relaxation) are allocated to the main direction
% (basically all DoFs are used to have the angular resolution
% however, if simply using a larger epsilon to relax more, the problem will
% only has a solution with smaller l1 norm, which is putting less energy
