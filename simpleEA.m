function [bestSoFarFit ,bestSoFarSolution ...
    ]=simpleEA( ...  % name of your simple EA function
    fitFunc, ... % name of objective/fitness function
    T, ... % total number of evaluations
    input) % replace it by your input arguments

% Check the inputs
if isempty(fitFunc)
  warning(['Objective function not specified, ''' objFunc ''' used']);
  fitFunc = 'objFunc';
end
if ~ischar(fitFunc)
  error('Argument FITFUNC must be a string');
end
if isempty(T)
  warning(['Budget not specified. 1000000 used']);
  T = '1000000';
end
eval(sprintf('objective=@%s;',fitFunc));
%% TODO
% Initialise variables
nbGen = 0; % generation counter
nbEval = 0; % evaluation counter
bestSoFarFit = 0; % best-so-far fitness value
bestSoFarSolution = NaN; % best-so-far solution
%recorders
fitness_gen=[]; % record the best fitness so far
solution_gen=[];% record the best phenotype of each generation
fitness_pop=[];% record the best fitness in current population 
%% Below starting your code
population_size = 4;
UPPER_BOUND = 31;
LOWER_BOUND = 0;
population = randi([LOWER_BOUND, UPPER_BOUND], 1, population_size);
population_gene = dec2bin(population, 5);
% Initialise a population
%% TODO
fitness = objFunc(population);

fitness_pop = [fitness_pop, max(fitness)];

for i = 1 : population_size
    if fitness(i) > bestSoFarFit
        bestSoFarFit = fitness(i);
        bestSoFarSolution = population(i);
    end
end

nbGen = nbGen + 1;
nbEval = nbEval + population_size;
fitness_gen = [fitness_gen, bestSoFarFit];
solution_gen = [solution_gen, bestSoFarSolution];


% bestSoFarFit = max([fitness, bestSoFarFit]);
% bestSoFarSolution = population;
% fitness_gen = [fitness_gen, bestSoFarFit];
% solution_gen = [solution_gen, population_gene];

% Evaluate the initial population
%% TODO

% Start the loop
while (nbEval<T) 
% Selection
prop = fitness./sum(fitness);

idx1 = population_size;
idx2 = population_size;
idx3 = population_size;
idx4 = population_size;
random_choice = rand;
tmp_choice = 0;
for i = 1:population_size - 1
    tmp_choice = tmp_choice + prop(i);
    if tmp_choice <= random_choice && tmp_choice + prop(i+1) > random_choice
        idx1 = i;
        break;
    end
end

fitness_tmp = fitness(idx1);
fitness(idx1) = 0;
prop = fitness./sum(fitness);

random_choice = rand;
tmp_choice = 0;
for i = 1:population_size-1
    tmp_choice = tmp_choice + prop(i);
    if tmp_choice <= random_choice && tmp_choice + prop(i+1) > random_choice
        idx2 = i;
        break;
    end
end

fitness(idx1) = fitness_tmp;

random_choice = rand;
tmp_choice = 0;
for i = 1:population_size-1
    tmp_choice = tmp_choice + prop(i);
    if tmp_choice <= random_choice && tmp_choice + prop(i+1) > random_choice
        idx3 = i;
        break;
    end
end

fitness_tmp = fitness(idx3);
fitness(idx3) = 0;
prop = fitness./sum(fitness);


random_choice = rand;
tmp_choice = 0;
for i = 1:population_size-1
    tmp_choice = tmp_choice + prop(i);
    if tmp_choice <= random_choice && tmp_choice + prop(i+1) > random_choice
        idx4 = i;
        break;
    end
end

fitness(idx3) = fitness_tmp;

random_cross_idx = randi(5);

new_gene_1 = horzcat(population_gene(idx1, 1:random_cross_idx), population_gene(idx2, random_cross_idx+1:5));
new_gene_2 = horzcat(population_gene(idx2, 1:random_cross_idx), population_gene(idx1, random_cross_idx+1:5));

random_cross_idx = randi(5);

new_gene_3 = horzcat(population_gene(idx3, 1:random_cross_idx), population_gene(idx4, random_cross_idx+1:5));
new_gene_4 = horzcat(population_gene(idx4, 1:random_cross_idx), population_gene(idx3, random_cross_idx+1:5));

%population_gene = [string(new_gene_1), string(new_gene_2), string(new_gene_3), string(new_gene_4)];
population_gene = [new_gene_1; new_gene_2; new_gene_3; new_gene_4];
% for i = 1: population_size
%     population_gene(i) = char(population_gene(i))
% end

% Reproduction (selection, crossver)

%% TODO

for i = 1: population_size
    if rand < 0.2
        rand_idx = randi(5);
        current_gene = population_gene(i,:);
        current_gene(rand_idx) ='1' - current_gene(rand_idx) + '0';
        population_gene(i,:) = current_gene;
    end
end



% Mutation
%% TODO
population = bin2dec(population_gene);

fitness = objFunc(population);

fitness_pop = [fitness_pop, max(fitness)];

for i = 1 : population_size
    if fitness(i) > bestSoFarFit
        bestSoFarFit = fitness(i);
        bestSoFarSolution = population(i);
    end
end


fitness_gen = [fitness_gen, bestSoFarFit];
solution_gen = [solution_gen, bestSoFarSolution];

nbGen = nbGen + 1;
nbEval = nbEval + population_size;
end
bestSoFarFit
bestSoFarSolution


figure,plot(1:nbGen,fitness_gen,'b') 
title('Fitness\_Gen')

figure,plot(1:nbGen,solution_gen,'b') 
title('Solution\_Gen')

figure,plot(1:nbGen,fitness_pop,'b') 
title('Fitness\_Pop')
%%
% fitness = [6,4,6,2,3,6]
% a = max(fitness(1:3))

