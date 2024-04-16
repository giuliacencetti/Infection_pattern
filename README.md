Code for the article https://arxiv.org/pdf/2309.10486.pdf.

The code language is Python 3, we make use of the following libraries: numpy, seaborn, matplotlib.pyplot, csv, networkx, pickle, os, sys, scipy.

We use Sociopatterns data, which are temporal, so first we need to aggregate them by using the notebooks:
- Generate_aggregated_graphs.ipynb (which, starting from data like Sociopatterns coming in a shape "i j t", generates static graphs)
- Generate_aggregated_hypergraphs.ipynb (which, starting from the same data, generates hypergraphs)
- Generate_aggregated_graphs_packformat.ipynb


**Simple model**:

- Simple_model_pack contains a simplified version of the code used in https://github.com/diegocontr/EpidemicSimulation to simulate simple contagion with R-calibrated betas.

- Simulate_simple_contagion_SIR_free_a.ipynb, which uses the functions in Simplagion_functions_SIR.py, simulates simple contagion and saves infection pattern in results.

- Simulate_simple_contagion_SIR_fixed_a_final.ipynb, which uses the functions in Simplagion_functions_SIR.py, simulates simple contagion and saves infection pattern in results only if the attack rate of the simulation is in the indicated range of attack rate (a_inf,a_sup).

- Simulate_simple_contagion_SIR_find_R0.ipynb, which uses the functions in Simple_functions_SIR_R0.py, simulates the contagion and compute R0.

- plot_similarity_matrix_SIR_simple_free_a.ipynb and plot_similarity_matrix_SIR_simple_fixed_a_final.ipynb compute and plot the cosine similarity between infection patterns at different values of beta.

- Receiver_spreader_index_simple_vs_simple.ipynb computes receiver and spreader indices and make the plots.


**Simplicial model**:

- Simulate_simplicial_contagion_SIR_free_a.ipynb, which uses the functions in Simplagion_functions_SIR.py, simulates simplicial contagion and saves infection pattern in results.

- Simulate_simplicial_contagion_SIR_fixed_a_final.ipynb, which uses the functions in Simplagion_functions_SIR.py, simulates simplicial contagion and saves infection pattern in results only if the attack rate of the simulation is in the indicated range of attack rate (a_inf,a_sup).

- Simulate_simplicial_contagion_SIR_find_R0.ipynb, which uses the functions in Simplagion_functions_SIR_R0.py, simulates the contagion and compute R0.

- plot_similarity_matrix_SIR_simplicial_free_a.ipynb and plot_similarity_matrix_SIR_simplicial_fixed_a_final.ipynb compute and plot the cosine similarity between infection patterns at different values of beta and beta_Delta.

- Receiver_spreader_index_simplicial_vs_simplicial_free_a.ipynb and Receiver_spreader_index_simplicial_vs_simplicial_fixed_a.ipynb compute receiver and spreader indices and make the plots.

**Threshold model**:

- Simulate Threshold_stochastic.ipynb, which uses the functions in Threshold_strength_responsibility_model_functions.py, simulates threshold contagion and saves infection pattern in results.

- plot_similarity_matrixSIR_threshold.ipynb computes and plots the cosine similarity between infection patterns at different values of theta.

- Receiver_spreader_index_threshold_vs_threshold.ipynb computes receiver and spreader indices and makes the plots.


**Models comparison**:

Notebooks that analyze and compare infection pattern across simple, simplicial and threshold model

