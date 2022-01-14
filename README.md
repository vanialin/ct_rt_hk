# ct_rt_hk
 
Codes for generating results in the manuscript "Incorporating temporal distribution of population-level viral load enables real-time estimation of COVID-19 transmissibility"

Data are not available to public due to data consent with relevant party but variable explanation has been written within each script (if necessary) to facilitate interpretation.

## main
Contains codes for generating main Fig. 1-2 and Supplementary Table 1-4.

* 1_merge_data - calculate daily Ct mean and skewness; merges Ct linelist and case count linelist. 
* 2_ct_for_bootstrap - use bootstrap to calculate the 95% CI for smoothed daily Ct (by GAM) and daily Ct skewness. 
* **3_plot_fig_1** - recreate Fig. 1 based on output data from these two scripts above. 
* 4_model_loglinear - daily Ct distribution regressed on incidence-based Rt to generate the model for Ct-predicted Rt. This script recreates Supplementary Tables 1~4. 
* **5_plot_fig_2** - recreate Fig. 2 based on output data from script 4 above. 

## supplementary
Contains codes for generating all Supplementary Figures as well as Supplementary Table 7.

### empirical
Here we demonstrated underlying mechanisms of Ct methods methods, conducted sensitivity analyses and further validated our model using empirical data in Hong Kong.
* s1_temporal_ct_delay - demonstrate temporal trends of daily Ct values and daily delay distribution of onset-to-sampling; reproduce Fig. S1.
* s2_ct_backprojection - backproject Ct values at onset for symptomatic cases and compared distributions of Ct at onset vs. Ct at sampling over time; reproduce Fig. S2-4.
* s3_model_sens_age - sensitivity analysis on temporal age distribution of confirmed cases on our estimates; reproduce Fig. S5.
* s4_model_sens_symp - sensitivity analysis to evaluate the impact of asymtomatic cases at detection on our estimates; reproduced Fig. S6.
* s5_model_reverse - use earlier stage in wave 4 to train the model; reproduce Fig. S7.
* s6_model_k-fold_validation - 10-fold cross validation on our model.
* s7_model_training_selection - evaluation of the performance of our model under various length and starts of training periods; reproduced Fig. S8.

### simulation
Here we conducted simulation recovery to demonstrate the utility of Ct-based method under different scenarios of case detection under symptom-based surveillance.
* sim_1_all_linelists - simulate incidence curve and all linelists under the 4 scenarios.
* sim_2_get_results - reproduce Fig. S9-10.
* sim_3_run_100times - run 100 times to get the CI for spearman's correlation rho between incidence-based Rt and Ct-Rt under each simulated scenario; reproduce Supplementary Tables 7.
* sim_funcs/sim_source - contain sourced functions used for simulation.   
In "partab_seir_switch_model_hk.csv" we provided parameters used for the transmission mode.


[Operated under R version 4.1.2 (R Development Core Team, 2021)]
