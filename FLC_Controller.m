% Custom Genetic Algorithm for Fuzzy Logic Controller Optimization

% Initialize the Fuzzy Inference System (FIS)
fis = mamfis('Name', 'FLC_Controller');

% Define the input variable "Temperature" with placeholder membership functions
fis = addInput(fis, [0 40], 'Name', 'Temperature');
fis = addMF(fis, 'Temperature', 'trimf', [0 0 15], 'Name', 'Low');  
fis = addMF(fis, 'Temperature', 'trimf', [10 20 30], 'Name', 'Medium');
fis = addMF(fis, 'Temperature', 'trimf', [25 40 40], 'Name', 'High');

% Define the input variable "Lighting" with initial membership functions
fis = addInput(fis, [0 1], 'Name', 'Lighting');
fis = addMF(fis, 'Lighting', 'trimf', [0 0 0.5], 'Name', 'Dim');
fis = addMF(fis, 'Lighting', 'trimf', [0.25 0.5 0.75], 'Name', 'Normal');
fis = addMF(fis, 'Lighting', 'trimf', [0.5 1 1], 'Name', 'Bright');

% Define the output variable "HVAC" with initial membership functions
fis = addOutput(fis, [0 1], 'Name', 'HVAC');
fis = addMF(fis, 'HVAC', 'trimf', [0 0 0.5], 'Name', 'Low');
fis = addMF(fis, 'HVAC', 'trimf', [0.25 0.5 0.75], 'Name', 'Medium');
fis = addMF(fis, 'HVAC', 'trimf', [0.5 1 1], 'Name', 'High');

% Define fuzzy rules for the system
ruleList = [ ...
    "If Temperature is Low and Lighting is Dim then HVAC is Low"; ...
    "If Temperature is Low and Lighting is Normal then HVAC is Medium"; ...
    "If Temperature is Low and Lighting is Bright then HVAC is Medium"; ...
    "If Temperature is Medium and Lighting is Dim then HVAC is Medium"; ...
    "If Temperature is Medium and Lighting is Normal then HVAC is Medium"; ...
    "If Temperature is Medium and Lighting is Bright then HVAC is High"; ...
    "If Temperature is High and Lighting is Dim then HVAC is Medium"; ...
    "If Temperature is High and Lighting is Normal then HVAC is High"; ...
    "If Temperature is High and Lighting is Bright then HVAC is High"];
fis = addRule(fis, ruleList);

% --- Simulation before optimization ---
inputValues = [5 0.2; 15 0.5; 30 0.8; 35 0.9];  % Different input values for testing
numTests = size(inputValues, 1);
outputsBefore = zeros(numTests, 1);

% Evaluate the initial FLC (before optimization)
for i = 1:numTests
    outputsBefore(i) = evalfis(fis, inputValues(i, :));
end

% Parameters for the genetic algorithm
populationSize = 20;
numGenerations = 50;
mutationRate = 0.2;
crossoverRate = 0.7;
numParams = 9;  % We are optimizing the Temperature MFs (Low, Medium, High)

lowerBounds = [0 0 0 10 20 25 25 40 40];  % Lower bound for Temperature MFs
upperBounds = [5 10 15 15 25 35 30 40 40];  % Upper bound for Temperature MFs

% Generate initial population
population = lowerBounds + (upperBounds - lowerBounds) .* rand(populationSize, numParams);
fitness = zeros(populationSize, 1);

% Evaluate fitness of initial population
for i = 1:populationSize
    fitness(i) = evaluateFitness(population(i,:), fis);
end

% Genetic algorithm loop
for gen = 1:numGenerations
    disp(['Generation ' num2str(gen) ' - Best Fitness: ' num2str(min(fitness))]);
    
    % Selection using tournament selection
    selected = tournamentSelection(population, fitness, populationSize);
    
    % Crossover
    offspring = selected;
    for i = 1:2:populationSize
        if rand() < crossoverRate && i+1 <= populationSize
            [offspring(i,:), offspring(i+1,:)] = crossover(selected(i,:), selected(i+1,:));
        end
    end
    
    % Mutation
    for i = 1:populationSize
        if rand() < mutationRate
            offspring(i,:) = mutate(offspring(i,:), lowerBounds, upperBounds);
        end
    end
    
    % Evaluate fitness of offspring
    for i = 1:populationSize
        fitness(i) = evaluateFitness(offspring(i,:), fis);
    end
    
    % Replace the population with the offspring
    population = offspring;
end

% Get best solution
[~, bestIdx] = min(fitness);
bestParams = population(bestIdx, :);

% Display the best parameters and fitness
disp('Best Parameters:');
disp(bestParams);
disp(['Best Fitness: ' num2str(min(fitness))]);

% Update FIS with the best parameters
fis = setMFParams(fis, bestParams);

% --- Simulation after optimization ---
outputsAfter = zeros(numTests, 1);

% Evaluate the optimized FLC (after optimization)
for i = 1:numTests
    outputsAfter(i) = evalfis(fis, inputValues(i, :));
end

% Plot and compare results before and after optimization
figure;
subplot(2,1,1);
bar([outputsBefore outputsAfter]);
legend('Before Optimization', 'After Optimization');
xlabel('Test Cases');
ylabel('HVAC Output');
title('Comparison of FLC Output Before and After Optimization');

% Plot the optimized FIS membership functions
figure;
subplot(3,1,1);
plotmf(fis, 'input', 1);
title('Optimized Temperature Membership Functions');
subplot(3,1,2);
plotmf(fis, 'input', 2);
title('Lighting Membership Functions');
subplot(3,1,3);
plotmf(fis, 'output', 1);
title('HVAC Membership Functions');

% --- Functions used in the GA ---
% Fitness evaluation function
function cost = evaluateFitness(params, fis)
    % Update FIS with the new parameters
    fis = setMFParams(fis, params);
    
    % Test input and desired output
    testInput = [25 0.75];  % Test input for evaluation
    desiredOutput = 0.75;  % Example desired output
    
    % Get actual output from the FIS
    actualOutput = evalfis(fis, testInput);
    
    % Calculate error as cost
    cost = abs(desiredOutput - actualOutput);
end

% Function to set membership function parameters with sorted values
function fis = setMFParams(fis, params)
    % Sort the parameters to ensure they are valid for 'trimf'
    lowParams = sort(params(1:3));  % Ensure parameters for 'Low' MF are sorted
    mediumParams = sort(params(4:6));  % Ensure parameters for 'Medium' MF are sorted
    highParams = sort(params(7:9));  % Ensure parameters for 'High' MF are sorted
    
    % Update the 'Temperature' membership functions with sorted parameters
    fis.Inputs(1).MembershipFunctions(1).Parameters = lowParams;  % Low
    fis.Inputs(1).MembershipFunctions(2).Parameters = mediumParams;  % Medium
    fis.Inputs(1).MembershipFunctions(3).Parameters = highParams;  % High
end

% Tournament selection
function selected = tournamentSelection(population, fitness, populationSize)
    selected = zeros(size(population));
    for i = 1:populationSize
        candidates = randi(populationSize, [2, 1]);
        if fitness(candidates(1)) < fitness(candidates(2))
            selected(i,:) = population(candidates(1), :);
        else
            selected(i,:) = population(candidates(2), :);
        end
    end
end

% Crossover function
function [child1, child2] = crossover(parent1, parent2)
    alpha = rand();
    child1 = alpha * parent1 + (1 - alpha) * parent2;
    child2 = alpha * parent2 + (1 - alpha) * parent1;
end

% Mutation function
function mutant = mutate(individual, lb, ub)
    mutationPoint = randi(length(individual));
    mutant = individual;
    mutant(mutationPoint) = lb(mutationPoint) + (ub(mutationPoint) - lb(mutationPoint)) * rand();
end
