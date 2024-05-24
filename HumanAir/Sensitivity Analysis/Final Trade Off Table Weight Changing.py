import numpy as np
import matplotlib.pyplot as plt

# Define the trade off table
criteria = ["Ground Clearance", "Emission reduction", "Development Risk", "Cost", "Operability"]
initial_weights = np.array([0.1, 0.3, 0.2, 0.3, 0.1])
concepts = ['Flying Wing', 'Canard', 'Conventional']

# Define importance
initial_importance = np.array([
    [3, 1, 1, 2, 1],
    [3, 4, 4, 4, 3],
    [4, 3, 5, 5, 5]
])


# Perform sensitivity analysis
weight_ranges = np.linspace(0, 1, 22)  # From 0 to 1 in steps of 0.05

# To store results for each criteria being varied
sensitivity_results = {concept: {criteria[i]: [] for i in range(len(criteria))} for concept in concepts}


# Iterate through each criteria weight
for i in range(len(initial_weights)):
    for w in weight_ranges:
        weights = initial_weights.copy()
        weights[i] = w
        weights /= np.sum(weights)  # Normalize to ensure weights sum to 1
        final_scores = np.sum(initial_importance * weights, axis=1)
        for j, concept in enumerate(concepts):
            sensitivity_results[concept][criteria[i]].append(final_scores[j])

# Plot the sensitivity results
fig, axes = plt.subplots(len(initial_weights), 1, figsize=(10, 15))

for i, criterion in enumerate(criteria):
    ax = axes[i]
    for concept in concepts:
        ax.plot(weight_ranges, sensitivity_results[concept][criterion], label=concept)
    ax.set_xlabel('Weight Value')
    ax.set_ylabel('Final Score')
    ax.set_title(f'{criterion}')
    ax.legend(loc = 'lower right')
    ax.grid()

plt.tight_layout()
plt.savefig("Sensitivity Analysis of Weight.pdf")
plt.show()



# Perform sensitivity analysis for importance
importance_ranges = range(1, 6)  # From 1 to 5
sensitivity_results_importance = {concept: {criteria[i]: [] for i in range(len(criteria))} for concept in concepts}

for i in range(initial_importance.shape[1]):
    for imp in importance_ranges:
        importance = initial_importance.copy()
        importance[:, i] = imp
        final_scores = np.sum(importance * initial_weights, axis=1)
        for j, concept in enumerate(concepts):
            sensitivity_results_importance[concept][criteria[i]].append(final_scores[j])

# Plot the sensitivity results for importance
fig, axes = plt.subplots(len(criteria), 1, figsize=(10, 15))

for i, criterion in enumerate(criteria):
    ax = axes[i]
    for concept in concepts:
        ax.plot(importance_ranges, sensitivity_results_importance[concept][criterion], label=concept)
    ax.set_ylabel('Final Score')
    ax.set_xlabel('Importance Value')
    ax.set_title(f'{criterion}')
    ax.legend(loc = 'lower right')
    ax.grid(True)

plt.tight_layout()
plt.savefig("Sensitivity Analysis of Importance.pdf")
plt.show()