## Updating Regional R0 using data from importations and subsequent transmission
This is a repository that holds all of the information necessary for the methods described in this paper **Paper Link Needed**.


## Getting the data and code
Clone this repository to make a local copy on your machine. You will then need to obtain the raw data needed for the analysis (note, the importation data will be fake data, unless you explicitly received it from Texas DSHS or myself). Go to this [link](https://doi.org/10.18738/T8/HYZ53B), and download the compressed folders. Unzip them and copy them into your local repository.

## Running the code
1. Double click on the .Rproj file `rnot_updater.Rproj`
2. (optional) Open the `generate_scam_list.R` to get monte carlo samples for R0 estimation parameters. This script takes a long time (>24 hours) to complete, so I've provided the generated samples in the `data_produced/vector_suitability` folder if you don't want to run it yourself. If you want to test it out for yourself though, change the number of runs to ~1-10, and it runs relatively quickly.
3. Open `analyze_temperature_data.R`, and run it through. This compiles both the historic and actual county temperatures for Texas into a usable format for the analysis, and stores them in the data folder. It takes about 20 minutes to run on my machine.
4. Run `calc_monthly_rnots.R` to get R0 distributions for each county and month.
5. Run the `calc_dispersion_table.R` file to obtain the dispersion parameter to be used for the fitting (should be 0.12)
      - If printer minimum is different from 0.12, substitute it into the `cpp_fitting_fxns.R` file in the `find_rnot_ods()` function in place of 0.12, though this should not be necessary.
6. Run the posterior scaling factor estimation to get posteriors at every importation time point for the scaling factor.
    - Run the `create_alpha_like_job_file.R` to create job file for running posterior estimation for all scenarios considered from the paper (temperature used, reporting rate, and secondary transmission number in November). This script creates a file in the `launcher` folder where each line is a "job" that can be run using the command line to call an R script. There are two options to subsequently run all jobs.
    1. I've provided the .slurm script `INSERT NAME HERE`, which is the script used to call the jobs using the TACC high performance computers. running the command `sbatch XXXXXX` will submit the job to TACC for running, though you will need to tweak the locations for saving files in the following scripts YYYY. You may also need to tweak the parameters for running on your own nodes depending on if you use TACC or some other high performance computing environment.
    2. If you don't have High performance computing capability to schedule jobs, you could alternatively run each of these lines manually from your own terminal, or create a different R script that accomplishes the same task. You will also need to change the saving location in the following scripts YYYY.

7.  Run `calc_final_posterior_rnots.R`, which will run single MCMC for each parameter set, and save posterior distributions for all county R0s in a usable format for plotting. It is currently setup to run on high performance computer, but could be easily tweaked to run locally. This does take a long time though if done on a local machine.     

8. run `fake_mcmc_dat_generator.R` to get the fake data used for the first figure
9. run `cty_sec_trans.R` to get data for last figure in ms
10. run `sim_test_likelihood_fxn.R` to get supplemental figure showing the likelihood vs simulation results.
11. run `ms_fig_creator.R`, which can be used to make all of the figures from the manuscript


