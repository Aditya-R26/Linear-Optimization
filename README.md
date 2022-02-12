# Linear-Optimization
Using the Two-Phase Algorithm in Linear Optimization, solved Pig Farming - Feed Optimization problem and Production Planning - Parts Manufacturing Optimization problem.
Algorithm uses the Two-Phase Simplex method, where it first detects if the Linear programming problem is feasible or not.
Gaussian reduction was used to check if the problem at hand has full row rank. If not, algorithm detected and removed the
redundant constraints.
• Bland’s rule was implemented to stop the algorithm from cycling. Application:
a) Pig Farming – Feed Optimization
• Optimized the quantities of the available types of feed (corn, tankage, alfalfa) that should be given to each pig to meet certain nutritional (carbohydrates, protein, vitamins) requirements at a minimum cost.
b) Production Planning – Part Manufacturing Optimization
• Optimized the quantities of 3 types of parts that need to be manufactured across each of the 3 branch plants under constraints of labor and equipment capacity to maximize profit.
