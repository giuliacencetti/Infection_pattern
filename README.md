Code for the article https://arxiv.org/pdf/2309.10486.pdf.

The code language is Python 3, we make use of the following libraries: numpy, seaborn, matplotlib.pyplot, csv, networkx, pickle, os, sys, scipy.

We use Sociopatterns data, which are temporal, so first we need to aggregate them by using the notebooks:
- *Generate_aggregated_graphs.ipynb*, which, starting from data like Sociopatterns coming in a shape "i j t", generates a static weighted graph.
- *Generate_aggregated_hypergraphs.ipynb*, which, starting from the same data, generates a weighted hypergraph.
- *Generate_aggregated_graphs_packformat.ipynb*

We then simulate three contagion models (simple, simplicial and threshold) on the static weighted networks and hypergraphs. For each model we analyze the results and finally we compare the results obtained for the three different models.

# Simple model:

- *Simple_model_pack* contains a simplified version of the code used in https://github.com/diegocontr/EpidemicSimulation to simulate simple contagion with R-calibrated betas.

- *Simulate_simple_contagion_SIR_free_a.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simple contagion and saves infection pattern in results. Notice that the code for simple model and simplicial model is the same, with the difference that for simple model we always set $\beta_{\Delta} = 0$.

- *Simulate_simple_contagion_SIR_fixed_a_final.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simple contagion and saves infection pattern in results only if the attack rate of the simulation is in the indicated range of attack rate (a_inf,a_sup).

- *Simulate_simple_contagion_SIR_find_R0.ipynb*, which uses the functions in *Simple_functions_SIR_R0.py*, simulates the contagion and compute R0.

- *plot_similarity_matrix_SIR_simple_free_a.ipynb* and *plot_similarity_matrix_SIR_simple_fixed_a_final.ipynb* compute and plot the cosine similarity between infection patterns at different values of beta.

- *Receiver_spreader_index_simple_vs_simple.ipynb* computes receiver and spreader indices and make the plots.


# Simplicial model:

- *Simulate_simplicial_contagion_SIR_free_a.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simplicial contagion and saves infection pattern in results.

- *Simulate_simplicial_contagion_SIR_fixed_a_final.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simplicial contagion and saves infection pattern in results only if the attack rate of the simulation is in the indicated range of attack rate (a_inf,a_sup).

- *Simulate_simplicial_contagion_SIR_find_R0.ipynb*, which uses the functions in *Simplagion_functions_SIR_R0.py*, simulates the contagion and compute R0.

- *plot_similarity_matrix_SIR_simplicial_free_a.ipynb* and *plot_similarity_matrix_SIR_simplicial_fixed_a_final.ipynb* compute and plot the cosine similarity between infection patterns at different values of $\beta$ and $\beta_{Delta}$.

- *Receiver_spreader_index_simplicial_vs_simplicial_free_a.ipynb* and *Receiver_spreader_index_simplicial_vs_simplicial_fixed_a.ipynb* compute receiver and spreader indices and make the plots.

# Threshold model:

- *Simulate Threshold_stochastic.ipynb*, which uses the functions in *Threshold_strength_responsibility_model_functions.py*, simulates threshold contagion and saves infection pattern in results.

- *plot_similarity_matrixSIR_threshold.ipynb* computes and plots the cosine similarity between infection patterns at different values of theta.

- *Receiver_spreader_index_threshold_vs_threshold.ipynb* computes receiver and spreader indices and makes the plots.


# Models comparison:

Notebooks that analyze and compare infection pattern across simple, simplicial and threshold model.

- *centrality_vs_r_s_indices.ipynb* 

- *plot_similarity_matrix_simplicial_vs_simple.ipynb* computes and plots the cosine similarity between infection patterns in simple model at varying $\beta$ and in simplicial model at varying ($\beta$, $\beta_{\Delta}$).

- *plot_similarity_matrix_threshold_vs_simple.ipynb* computes and plots the cosine similarity between infection patterns in simple model at varying $\beta$ and threshold model at varying $\theta$.

- *Receiver_spreader_index_simple_vs_simplicial_free_a.ipynb* computes receiver and spreader indices for simple and simplicial models and makes the comparison plots.

- *Receiver_spreader_index_simple_vs_simplicial_fixed_a.ipynb* does the same but with fixed attack rate.

- *Receiver_spreader_index_simple_vs_threshold.ipynb* computes receiver and spreader indices for simple and threshold models and makes the comparison plots.

