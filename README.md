# CHEME5999
Independent study aimed at modeling chip-based glycosylation by bacterial oligosaccharyltransferase 
Program Instruction:
====================

The following instructions outline those described in CHEME5999_writeup_FKP.pdf   All cade is written in Matlab and was run using version R2018a. The governing ODEs are written to describe species concentration over time in a microfluidic device where enzymatic glycosylation of a target protein is performed.  

Our goal was to (A) develop a model for chip-based glycosylation  and (B) to explore different methods for estimating key parameters, including: 
   Km1 - michaelis mention constant for sugar nucleotide
   Km2 - michaelis mention constant for bare peptide, 
   kcat  -  turnover rate 
   n - cooperativity   

The final program  CHEME5999_FinalPset.m compares model predictions based on an original proposed parameter set  to simulations performed with (1) real parameters taken from literature, (2) parameters estimated using a simulated annealing (SA), and (3) parameters learned from “fake” experiments.   
Comparison of original parameter versus “real”(1) or “learned” parameters (2&3) 

Run CHEME5999_FinalPset.m . The script will generate 3 graphs comparing the original model to  model predictions using real parameters (fig 1),  SA parameters (fig 2), and experimental parameters (fig 3). 

1.  Real Parameters :  real parameters were taken from Gerber, S. et al. Mechanism of Bacterial Oligosaccharyltransferase: In vitro quanitification of sequon binding and caralysis. J. Biol. Chem. 2013, 288, 8849 – 8861. 

2.  Simulated Annealing Parameters: CHEME5999_PS2_SA.m  is a learning algorithm that uses simulated annealing to estimate Km1, Km2, kcat, and n

3. Experimental Design Parameters. See  “Experimental Design Algorithm”  section and associated folder “ CHEME5999_ED_Algorithm”.  This set of files uses fake experiments to back out estimates for Km1, Km2, kcat, and n.   


--------------------------------------------------------------------------------------------------------------
Experimental Design Algorithm
This section describes how to use the following compilation of matlab files to estimate values for Km1, Km2, kcat, and n.  The related files include:
        1. EXP_Data_x1.m             11. Sensitivity Analysis         21. Estimation_Param.m
        2. EXP_Data_x2.m             12. DataFile.m                   22. Estimation_fmincon.m 
        3. EXP_Data_x3.m	           13. Exp_Data.m		                23. P_solution2.txt
        4. EXP_Data__Eo.m	           14. Call_ODE.m                   24. Analysis_Param.m
        5. EXP_Data_x12.m            15. Analytical_SSM.m             25. Analysis_func.m 
        6. Data_EXP1.txt	           16. senfunction.m                26. DATA_EXP_3.txt 
        7. Data_EXP2.txt             17. Jacobian.m
        8. Data_EXP3.txt             18. P_matrix.m 
        9. Data_EXP4.txt             19. Normalize_SSM.m
       10. Data_EXP5.txt             20. Identifiability.m 

I.	Generate Fake Experimental Data
        a.	EXP_Data_x1/x2/x3/Eo/x12.m ....  Data_EXP1/2/3/4/5.txt
            i. EXP_Data_x1/x2/x3/Eo/x12.m (script) :  different fake experiments
                 1. Each experiment alters initial concentration and guess for each state (i.e. species) variable.
                 2. Saves output in txt file 
           ii. Data_EXP1/2/3/4/5.txt: data output 
                 1. Contains concentration (of each species) vs. time data for each experiment. 

II.	Perform sensitivity analysis. Run Sensitivity_Analysis.m (script) and the ‘3” when prompted for experiment number (EXP_NM). I only used experiment 3 to estimate parameters.  
        a.	DataFile.m (fxn) 
             i.	For each experiment, specifies both the initial guesses for each state variable and the enzyme concentration in the device
            ii. Specifies initial “guess” parameter set (based on model from pset 1)
            iii. Defines range of guesses for optimal parameter sets
        b.	Exp_Data.m (fxn) 
              i.Loads experimental data from text files generated in STEP I, and passes them to main script for sensitivity analysis (Sensitivty_Analysis.m). 
        c. 	Call_ODE.m (fxn) 
              i.	Defines model equations and conditions for each experiment (exact same conditions specified in EXP_DATA_x1/x2/…m scripts. 
        d.	Analytical_SSM.m (fxn) 
              i.	Computes sensitivity matrix  ( “z” matrix) using sensitivity coefficient model equations (senfunction.m, see below), which depends on the Jacobian matrix (Jacobian.m, see below) and  the Parameter Matrix ( P_matrix.m, see below). 
              ii. Jacobian.m (fxn)
                    1.	Defines  Jacobian matrix derived from model equations. 
                    2.	Specifies variations in parameter values between experiments, specifically  the Enzyme concentration. 
             iii. P_matrix.m
                    1.	Defines Parameter Matrix as derived from model equations 
                    2.	Specifies variation in the parameter values (enzyme concentration) and in the initial values of state variables between experiments.   
              iv. senfunction.m 
                    1.	Function  defines model equations for computing sensitivity coefficient (z) matrix
               v. Normalize_SSM.m (fxn) 
                    1.	Normalize sensitivity matrix (SSM) ouput by  senfunction.m 
        e.	Identifiability.m (fxn)
               i.	 Determines which parameters are identifiable based on SSM computed in Analytical_SSM.m 
              ii.	 Follows McAuley procedure (doi:10.1081 PRE-120024426) 

III. Parameter estimation. Run Estimation _Param.m (script). Note, that estimations  are only performed using Experiment 3. 
        a.	DataFile.m (fxn) 
              i.	See part II for description 
        b.	Exp_Data.m (fxn)
              i.	See part II for description 
        c.	Call_ODE,m  (fxn) 
              i.	See part II for description 
        d.	Estimation_fmincon.m (fxn)
              i.	Uses fmincon as learning algorithm to estimate optimal parameter values.  
        e.	Writes and saves estimated parameter sets for all iterations to P_solution2.txt 
        f. I used the last row of numbers in the P_solution2.txt file as my optimized parameters. These numbers were used to generate figure 3 in the CHEME5999_FinalPset.m 
              i.  NOTE, the parameters are number as follows:   P(1) = KM1 , P(2): = KM2, P(3) = kcat, P(4) = n 

IV.	Analyze estimated parameters. Run Analysis_Param.m (script) 
        a.	Analysis_func.m (fxn)
              i.	Reads estimate parameters from P_Solution2.txt and passes  them to Analysis _Param.m to recreate trajectories with estimated parameters.  
        b.	Exp_Data.m (fxn)
              i.	See part II for description 
        c.	Call_ODE.m (fxn)
              i.	See part II for description
        d.	Write and saves graph of data points to DATA_EXP_3.txt 
              i.	The data points include mean, lower bound, upper bound,  and normalized experimental data versus time.    
