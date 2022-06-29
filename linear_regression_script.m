clc,clear all, close all


%% Load data
data = load('ex1data1.txt');
X = data(:, 1); y = data(:, 2);
m = length(y); % number of training examples

% Plot Data
% Note: You have to complete the code in plotData.m
plot(X, y, 'rx', 'MarkerSize', 10); % Plot the data
ylabel('Profit in $10,000s'); % Set the y-axis label
xlabel('Population of City in 10,000s'); % Set the x-axis label


X = [ones(m,1), X];

%% SGD 
teta = zeros(2,1);
iteration = 1500;
lr = 0.01;
J = zeros(iteration,1);

for i=1:iteration
    
    teta0 = teta(1) - lr*(1/m)*sum((X*teta - y));
    teta1 = teta(2) - lr*(1/m)*sum((X*teta - y)'*X(:,2));
    J(i,1) = (1/(2*m))*sum(((X*[teta0; teta1]) - y).^2);
    teta = [teta0;teta1];
    
end
 %%%% -3.6303
 %%%% 1.1664
 
 
 %%
 %% ============= Part 4: Visualizing J(theta_0, theta_1) =============
fprintf('Visualizing J(theta_0, theta_1) ...\n')

% Grid over which we will calculate J
theta0_vals = linspace(-10, 10, 100);
theta1_vals = linspace(-1, 4, 100);

% initialize J_vals to a matrix of 0's
J_vals = zeros(length(theta0_vals), length(theta1_vals));

% Fill out J_vals
for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];
	  J_vals(i,j) = (1/(2*m))*sum(((X*t) - y).^2);
    end
end


% Because of the way meshgrids work in the surf command, we need to
% transpose J_vals before calling surf, or else the axes will be flipped
J_vals = J_vals';
% Surface plot
figure;
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');

% Contour plot
figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 3, 20))
xlabel('\theta_0'); ylabel('\theta_1');
hold on;
plot(teta(1), teta(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
