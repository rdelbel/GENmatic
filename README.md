#GENmatic

##Fast and Automated genetic reporting and analysis 

GENmatic provides a consistent R interface to do genetic QC and Analysis. The following Genetic anlysis are currently possible
* coxph with covariates and stratification, additive, dominant, recessive model
* linear with covariates, additive model
* logistic with covariates, additive model
coxph is limited to candidate gene studies on windows.

###GENfit
The main analysis function is GENfit, which completes the procedures for an entire genentic association study. There are many arugments to the function, but some are optional. See the function documentation in R for details. GENfit calls many functions defined in GENmatic. These functions could also be used independently if you wish.
* make_covar_pheno makes covar and pheno file for plink not(windows and coxph)
* makeRplink makes R function to throw into plink to do cox model not(windows) and coxph
* plinkfit fits models in plink not(windows and coxph)
* GENwincox fits cox models (windows and coxph)
  * plinktoR converts plink data files to R dataframe
* snpinfo gets supplimentary info for each snp from plink bfile
* manhattan plots manhattan plot of results
* qq plots qq plot of results
* topsnps gets information about top snps from results

###How to use on on server
* Login to scinet in putty. login.scinet.utoronto.ca
* Put in username/password
* Login to development node $ssh gpcxx where xx is a number between 01 and 04
* One Time Setup:
  * Create a new folder GENmatic in your scratch drive $mkdir $SCRATCH/GENmatic
  * Transfer runGENmatic.sh to $SCRATCH/GENmatic (eg using scp or filezilla)
  * With filezilla select port 22 and hit quick connect
  * $chmod 711 $SCRATCH/GENmatic/runGENmatic.sh
* Make a new directory to hold all files from the analysis eg $mkdir $SCRATCH/MyNewGwas
* Transfer .bim/.bed/.fam and .csv corresponding to the clinical data into the new folder
* Open GENmatic.R and do not tuch the top code (unless you know what you are doing)
* Change pd= to be the location of your plink file on scinet 
* Write all the calls to GENfit you wish run
* Upload this file to the folder you creaded (eg $SCRATCH/MyNewGWAS)
* To run the code on the server $$SCRATCH/GENmatic/runGENmatic.sh <foldername> <plink filename> <data filename>
* In this case $$SCRATCH/GENmatic/runGENmatic.sh MyNewGwas thinned clin.csv
* $ $SCRATCH/runGENmatic.sh foldername GENmatic.R
Note that runGENmatic.sh will only let your code run for at most 2 hours. If you need more then you will need to edit runGENmatic.sh by changing walltime (up to 24h) and the queue. You can achieve this by changing walltime to your desired time and changing debug to batch. The longer you set the walltime the longer you will have to wait in the queue for available resources to be found.

###Future Work

There is still plenty of work that could be done.
* Improve the manhattaon plot
* Improve the qq plot
* Allow donominant, recessive linear/logistic models
* Interface with more advanced genetic software such as phase2

