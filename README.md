Code for the article https://arxiv.org/pdf/2309.10486.pdf.

The code language is Python 3, we make use of the following libraries: numpy, seaborn, matplotlib.pyplot, csv, networkx, pickle, json, os, sys, scipy.

We use Sociopatterns data, which are temporal, so first we need to aggregate them by using the notebooks:
- *Generate_aggregated_graphs.ipynb*, which, starting from data like Sociopatterns coming in a shape "i j t", generates a static weighted graph.
- *Generate_aggregated_hypergraphs.ipynb*, which, starting from the same data, generates a weighted hypergraph.

We then simulate three contagion models (simple, simplicial and threshold) on the static weighted networks and hypergraphs. For each model we analyze the results and finally we compare the results obtained for the three different models.

# Simple model:
The code for simple contagion is in the folders "Simple_model" and "Simple_model_pack".

## Simple_model_pack:

contains a simplified version of the code used in https://github.com/diegocontr/EpidemicSimulation to simulate simple contagion with R-calibrated betas.

To generate simulation data, execute the following scripts:

\```bash
python 01_sistematic_simulation_altmodel.py --disease {MODEL} --tauE {tauE} --tauI {tauI} --pvar {MARKOVIAN} --Nruns 1000 --R0fit {R0} --data {DATASET} --folder SimResults/{DATASET} --label {NAME}_ --modelname {NAME}
\```

Here's what each parameter represents:

- `MODEL`: Choose from SIR, SEIR, or COVID.
- `tauE`: Either 1 or 4, corresponding to 1/mu_E for the SEIR model.
- `tauI`: Either 1 or 4, corresponding to 1/mu_I for the SEIR and SIR models.
- `MARKOVIAN`: 1 for a Markovian process, or 0.25 for the non-Markovian process used in the paper.
- `R0`: Choose from 2, 2.5, 3, 3.5, or 4.0.
- `DATASET`: Select from office, hospital, cprepa, conf, or school.
- `NAME`: Use one of the following to ensure coherence with the parameters:
  - SIR: for the SIR model
  - COVID: for the COVID model
  - SEIRe1: for Markovian SEIR with tauI=1 and tauE=1
  - SEIRe4: for Markovian SEIR with tauI=1 and tauE=4
  - SEIRi4: for Markovian SEIR with tauI=4 and tauE=1
  - SEIRe1v025: for Non-Markovian SEIR with tauI=1 and tauE=1
  - SEIRe4v025: for Non-Markovian SEIR with tauI=1 and tauE=4
  - SEIRi4v025: for Non-Markovian SEIR with tauI=4 and tauE=1

For the school dataset, use the *01_sistematic_simulation_altmodel-school.py* script specifically for the COVID model.

To obtain the files needed to make the figures, use the script *02_extract_matrices.py*


## Simple_model:
Contains the following notebooks:

- *Simulate_simple_contagion_SIR_free_a.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simple contagion and saves infection pattern in results. Notice that the code for simple model and simplicial model is the same, with the difference that for simple model we always set $\beta_{\Delta} = 0$.

- *Simulate_simple_contagion_SIR_fixed_a_final.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simple contagion and saves infection pattern in results only if the attack rate of the simulation is in the indicated range of attack rate (a_inf,a_sup).

- *Simulate_simple_contagion_SIR_find_R0.ipynb*, which uses the functions in *Simple_functions_SIR_R0.py*, simulates the contagion and compute R0.

- *plot_similarity_matrix_SIR_simple_free_a.ipynb* and *plot_similarity_matrix_SIR_simple_fixed_a_final.ipynb* compute and plot the cosine similarity between infection patterns at different values of beta.

- *Receiver_spreader_index_simple_vs_simple.ipynb* computes receiver and spreader indices and make the plots.


# Simplicial model:
In the folder "Simplicial_model/code" there is the code for simplicial contagion, which uses the following notebooks:

- *Simulate_simplicial_contagion_SIR_free_a.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simplicial contagion and saves infection pattern in results.

- *Simulate_simplicial_contagion_SIR_fixed_a_final.ipynb*, which uses the functions in *Simplagion_functions_SIR.py*, simulates simplicial contagion and saves infection pattern in results only if the attack rate of the simulation is in the indicated range of attack rate (a_inf,a_sup).

- *Simulate_simplicial_contagion_SIR_find_R0.ipynb*, which uses the functions in *Simplagion_functions_SIR_R0.py*, simulates the contagion and compute R0.

- *plot_similarity_matrix_SIR_simplicial_free_a.ipynb* and *plot_similarity_matrix_SIR_simplicial_fixed_a_final.ipynb* compute and plot the cosine similarity between infection patterns at different values of $\beta$ and $\beta_{Delta}$.

- *Receiver_spreader_index_simplicial_vs_simplicial_free_a.ipynb* and *Receiver_spreader_index_simplicial_vs_simplicial_fixed_a.ipynb* compute receiver and spreader indices and make the plots.

# Threshold model:
In the folder "Threshold_model/code" there is the code for threshold contagion, which uses the following notebooks:

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


# How to reproduce the figures of the paper:
For all the figures we have simulated the processes on the primary school dataset of Sociopatterns.

- Fig. 3(a): Diego: explain how to generate the infection pattern matrix Cij.
Then the notebook *Plot_similarity_fixed_R0.ipynb* in "Simple_contagion/code" compares the generated infection patterns and generates the figure.

- Fig. 3(b): We first need to choose some values of $\beta$ (and hence of $R_0$) to simulate the simple contagion. For the figure in the paper the values of $\beta$ are [0.12,0.15,0.18,0.21,0.24,0.27,0.3,0.33,0.36], corresponding to the values of $R_0$ [1.40,1.65,1.97,2.11,2.31,2.48,2.68,2.84,3.01] (but different values can be found using *Simulate_simple_contagion_SIR_find_R0.ipynb*). These values can be put in *Simulate_simple_contagion_SIR_free_a.ipynb* in folder "Simple_model/code" as "beta_range" to simulate the process. The infection patterns will be stored in a folder "results".  The notebook *plot_similarity_matrix_SIR_simplicial_free_a.ipynb* takes these results and compare them to generate the figure.

- Fig.3(c-e): Similarly to the procedure for fig. 3(b), we need to simulate the process with *Simulate_simple_contagion_SIR_fixed_a_final.ipynb* and to make the plots with *plot_similarity_matrix_SIR_simple_fixed_a_final.ipynb*. In both the notebooks the desired range of attack rate should be indicated setting minimum (a_inf) and maximum (a_sup).

- Fig.4: ..Diego: explain how to generate the results.
Then the notebook *Fig_similarity_time.ipynb* in "Simple_model/code" will generate the figure.

- Fig.5(a,b): The notebook *Simulate_simplicial_contagion_SIR_free_a.ipynb* in "Simplicial_model/code" simulates the simplicial contagion for the set of values of $\beta$ and $\beta_{\Delta}$ that appears in the figure (but can be changed setting beta_betaT_range). The infection pattern are stored in the folder "results" and then used by the notebook *plot_similarity_matrix_SIR_simplicial_free_a.ipynb* to generate the figures.

- Fig. 5(c,d): Similarly to figs. 5(a,b), we need to simulate the process with *Simulate_simple_contagion_SIR_fixed_a_final.ipynb* and to make the figures with *plot_similarity_matrix_SIR_simplicial_fixed_a_final.ipynb*.  In both the notebooks the desired range of attack rate should be indicated setting minimum (a_inf) and maximum (a_sup).

- Fig. 5(e): The figure is generated by *plot_similarity_matrix_simplicial_vs_simple.ipynb* in "Models_comparison", which takes in input the infection patterns of simple and simplicial contagion generated with *Simulate_simple_contagion_SIR_free_a.ipynb* and *Simulate_simplicial_contagion_SIR_free_a.ipynb* respectively.

- Fig. 6(a): The notebook *Simulate Threshold_stochastic.ipynb* in "Threshold_model/code" simulates the process for the set of values of $\theta$ that appears in the figure, then *plot_similarity_matrixSIR_threshold.ipynb* generates the figure.

- Fig. 6(b): The figure is generated by *plot_similarity_matrix_threshold_vs_simple.ipynb* in "Models_comparison", which takes in input the infection patterns of simple and threshold contagion generated with *Simulate_simple_contagion_SIR_free_a.ipynb* and *Simulate Threshold_stochastic.ipynb* respectively.
