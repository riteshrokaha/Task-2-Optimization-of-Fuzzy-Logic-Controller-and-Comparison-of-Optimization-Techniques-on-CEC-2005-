% Comparison of Optimization Techniques on CEC'2005 Benchmark Functions
% Sphere Function and Rastrigin Function in 2D and 10D spaces
% Using Genetic Algorithm, Particle Swarm Optimization, and Simulated Annealing

clc;
clear;

% Define the dimensionalities
dim2D = 2;
dim10D = 10;
numRuns = 15;

% Store results for comparison
resultsGA_2D_Sphere = zeros(numRuns, 1);
resultsGA_10D_Sphere = zeros(numRuns, 1);
resultsPSO_2D_Sphere = zeros(numRuns, 1);
resultsPSO_10D_Sphere = zeros(numRuns, 1);
resultsSA_2D_Sphere = zeros(numRuns, 1);
resultsSA_10D_Sphere = zeros(numRuns, 1);

resultsGA_2D_Rastrigin = zeros(numRuns, 1);
resultsGA_10D_Rastrigin = zeros(numRuns, 1);
resultsPSO_2D_Rastrigin = zeros(numRuns, 1);
resultsPSO_10D_Rastrigin = zeros(numRuns, 1);
resultsSA_2D_Rastrigin = zeros(numRuns, 1);
resultsSA_10D_Rastrigin = zeros(numRuns, 1);

% Set optimization options
optsGA = optimoptions('ga', 'MaxGenerations', 100, 'PopulationSize', 50, 'Display', 'off');
optsPSO = optimoptions('particleswarm', 'MaxIterations', 100, 'SwarmSize', 50, 'Display', 'off');
optsSA = optimoptions('simulannealbnd', 'MaxIterations', 100, 'Display', 'off');

% ----- Sphere Function -----
disp('Running optimization on Sphere function...');
for i = 1:numRuns
    % Genetic Algorithm on Sphere (2D)
    [~, fval] = ga(@sphereFunction, dim2D, [], [], [], [], -5 * ones(1, dim2D), 5 * ones(1, dim2D), [], optsGA);
    resultsGA_2D_Sphere(i) = fval;

    % Genetic Algorithm on Sphere (10D)
    [~, fval] = ga(@sphereFunction, dim10D, [], [], [], [], -5 * ones(1, dim10D), 5 * ones(1, dim10D), [], optsGA);
    resultsGA_10D_Sphere(i) = fval;

    % Particle Swarm Optimization on Sphere (2D)
    [~, fval] = particleswarm(@sphereFunction, dim2D, -5 * ones(1, dim2D), 5 * ones(1, dim2D), optsPSO);
    resultsPSO_2D_Sphere(i) = fval;

    % Particle Swarm Optimization on Sphere (10D)
    [~, fval] = particleswarm(@sphereFunction, dim10D, -5 * ones(1, dim10D), 5 * ones(1, dim10D), optsPSO);
    resultsPSO_10D_Sphere(i) = fval;

    % Simulated Annealing on Sphere (2D)
    [~, fval] = simulannealbnd(@sphereFunction, rand(1, dim2D), -5 * ones(1, dim2D), 5 * ones(1, dim2D), optsSA);
    resultsSA_2D_Sphere(i) = fval;

    % Simulated Annealing on Sphere (10D)
    [~, fval] = simulannealbnd(@sphereFunction, rand(1, dim10D), -5 * ones(1, dim10D), 5 * ones(1, dim10D), optsSA);
    resultsSA_10D_Sphere(i) = fval;
end

% ----- Rastrigin Function -----
disp('Running optimization on Rastrigin function...');
for i = 1:numRuns
    % Genetic Algorithm on Rastrigin (2D)
    [~, fval] = ga(@rastriginFunction, dim2D, [], [], [], [], -5 * ones(1, dim2D), 5 * ones(1, dim2D), [], optsGA);
    resultsGA_2D_Rastrigin(i) = fval;

    % Genetic Algorithm on Rastrigin (10D)
    [~, fval] = ga(@rastriginFunction, dim10D, [], [], [], [], -5 * ones(1, dim10D), 5 * ones(1, dim10D), [], optsGA);
    resultsGA_10D_Rastrigin(i) = fval;

    % Particle Swarm Optimization on Rastrigin (2D)
    [~, fval] = particleswarm(@rastriginFunction, dim2D, -5 * ones(1, dim2D), 5 * ones(1, dim2D), optsPSO);
    resultsPSO_2D_Rastrigin(i) = fval;

    % Particle Swarm Optimization on Rastrigin (10D)
    [~, fval] = particleswarm(@rastriginFunction, dim10D, -5 * ones(1, dim10D), 5 * ones(1, dim10D), optsPSO);
    resultsPSO_10D_Rastrigin(i) = fval;

    % Simulated Annealing on Rastrigin (2D)
    [~, fval] = simulannealbnd(@rastriginFunction, rand(1, dim2D), -5 * ones(1, dim2D), 5 * ones(1, dim2D), optsSA);
    resultsSA_2D_Rastrigin(i) = fval;

    % Simulated Annealing on Rastrigin (10D)
    [~, fval] = simulannealbnd(@rastriginFunction, rand(1, dim10D), -5 * ones(1, dim10D), 5 * ones(1, dim10D), optsSA);
    resultsSA_10D_Rastrigin(i) = fval;
end

% ----- Displaying Results -----

% Compute mean results
meanGA_Sphere = [mean(resultsGA_2D_Sphere), mean(resultsGA_10D_Sphere)];
meanPSO_Sphere = [mean(resultsPSO_2D_Sphere), mean(resultsPSO_10D_Sphere)];
meanSA_Sphere = [mean(resultsSA_2D_Sphere), mean(resultsSA_10D_Sphere)];

meanGA_Rastrigin = [mean(resultsGA_2D_Rastrigin), mean(resultsGA_10D_Rastrigin)];
meanPSO_Rastrigin = [mean(resultsPSO_2D_Rastrigin), mean(resultsPSO_10D_Rastrigin)];
meanSA_Rastrigin = [mean(resultsSA_2D_Rastrigin), mean(resultsSA_10D_Rastrigin)];

% Print results to the command window
disp('--- Sphere Function Results ---');
disp(['GA (2D): ' num2str(meanGA_Sphere(1)) ', GA (10D): ' num2str(meanGA_Sphere(2))]);
disp(['PSO (2D): ' num2str(meanPSO_Sphere(1)) ', PSO (10D): ' num2str(meanPSO_Sphere(2))]);
disp(['SA (2D): ' num2str(meanSA_Sphere(1)) ', SA (10D): ' num2str(meanSA_Sphere(2))]);

disp('--- Rastrigin Function Results ---');
disp(['GA (2D): ' num2str(meanGA_Rastrigin(1)) ', GA (10D): ' num2str(meanGA_Rastrigin(2))]);
disp(['PSO (2D): ' num2str(meanPSO_Rastrigin(1)) ', PSO (10D): ' num2str(meanPSO_Rastrigin(2))]);
disp(['SA (2D): ' num2str(meanSA_Rastrigin(1)) ', SA (10D): ' num2str(meanSA_Rastrigin(2))]);

% ----- Plotting Results -----
figure;
bar([meanGA_Sphere; meanPSO_Sphere; meanSA_Sphere]);
legend('2D', '10D');
xlabel('Optimization Technique');
ylabel('Mean Fitness Value');
xticklabels({'GA', 'PSO', 'SA'});
title('Comparison of Optimization Techniques on Sphere Function');

figure;
bar([meanGA_Rastrigin; meanPSO_Rastrigin; meanSA_Rastrigin]);
legend('2D', '10D');
xlabel('Optimization Technique');
ylabel('Mean Fitness Value');
xticklabels({'GA', 'PSO', 'SA'});
title('Comparison of Optimization Techniques on Rastrigin Function');


% ---- Benchmark Functions ----
% Sphere function (F1)
function y = sphereFunction(x)
    y = sum(x.^2);
end

% Rastrigin function (F7)
function y = rastriginFunction(x)
    n = length(x);
    A = 10;
    y = A * n + sum(x.^2 - A * cos(2 * pi * x));
end
