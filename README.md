# OPA-design
 This repository is for my thesis related works

## Summary
* The "src" contains the code for this project. 
* There is a function calculating and plotting the array factor named "AF"
* The script "reproducing_plot.m" plots three figures in "Thinned Arrays Using Genetic Algorithms" published on IEEE Transactions on Antennas and Propagation, by Randy L. Haupt. 
* By comparing the plots in paper and the plots made by "reproducing_plot.m", the function "AF" is verified to be correct. 

# Algorithms
## 1. Genetic algorithm

### 1.1 General introducrion

The Genetic algorithm contains the following steps
1. Set/encode the parameters into gene (The activation weight 0/1 in our case)
2. Randomly generate M number of genes. 
3. Fitness computation (The max side lobe level in our case)
4. Natural selection (Based on max side lobe level. the smaller, the better)
5. Paring and mating (Crossover)
6. Mutation (Important for avoiding local minimum)
```
START
Generate the initial population
REPEAT
    1. Compute fitness
    2. Selection
    3. Crossover
    4. Mutation 
UNTIL population has converged
STOP
```
### 1.2 Reproducing the result

The GA converged to the result from the paper conditionally. 

The condition is including the result from the reference paper in the initialization. Since the algorithm starts from random initializations and performs the evolution stochastically, it is very difficult to generate the exact same result with the reference paper. "Cheating" by including the exact result in the initialization, the algorithm successfully converged to it. The convergence means the result from the reference is at least a local minimum. 
(The function "disp_fr" can show the filling rate of each gene)

Since the mutation helps to jump out of the local minimum, I decided to use only mutation to generate new genes. (GA_intentively_mutation.m)

```
START
Generate the initial population
REPEAT
    1. Compute fitness
    2. Selection
    3. Mutation
UNTIL population has converged
STOP
```

By using this more aggressive mutation strategy, I found a solution with a 73 filling rate which is better than the one from the reference (stored in "wopt.m"). This might prove that the result from the reference is indeed a local minimum (At least there is a better solution). The script "Compare2.m" compares two results. 
However, this algorithm failed to converge. 

Because the paper found a local minimum and the algorithm converges to that local minimum with certain initialization, I might conclude that the GA has been correctly implemented. 
### 1.3 Discussion

The GA is intrinsically a combinatorial optimization, where the algorithm searches for an optimal solution within a finite set. 

GA starts with random initializations, the process of evolution is stochastic as well. Compared with the gradient descent method, where the step of each iteration is analytically calculated based on the gradient of the objective function, the GA is less interpretable. For the gradient descent method, the "route" of convergence is the same if the initial point is the same. The GA, however, has two steps introducing randomness, crossover and mutation, making reproducing more difficult. 

Therefore, I would like to conclude that the unconditionally exact reproduction of the result of the paper is hardly possible. 