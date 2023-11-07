Title: Optimization program and numerical data of simultaneous modeling of in vivo and in vitro effects of nondepolarizing neuromuscular blocking drugs

Author: Hikaru Hoshino and Eiko Furutani

Description: This data provides python source codes of optimization program for simultaneous modeling of in vivo and in vitro effects of Nondepolarizing
Neuromuscular Blocking Drugs (NDNBs) and numerical data of the optimization results of parameter estimation for three NDNBs of Cisatracurium, Vecuronium, and Rocuronium. 

License: CC BY 4.0 

---------
Contents: 

minimization.py: Main program of parameter estimation.
                 The result of parameter estimation is stored in "minimized_parameters_xx.py" as a python script,
                 and log information is stored in "minimization_xx.log". 
                 (xx is a model type)

minimization_kdissD.py: Program for parameter estimation with different kdissD1 and kdissD2. 

fig_effect_curves.py: Program for plotting concentration v.s. effect curves (Figure 3)

fig_parameter_sweep.py: Program for calculating pharmacologic parameters (EC50, IC50 etc) under various model parameters. 
                        The results are sotred in "fig_(model_type)/parameter_sweep_(invitro/invivo).csv". 
                        
fig_parameter_sweep_plot.py: Program for plotting the avobe results (Figure 5-7)

fig_time_course.py: Program for plotting the time course of activation AChRs (Figure 4)

modules: several utility functions and programes used by the above main programs.  
