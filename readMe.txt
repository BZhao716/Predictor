                        Installation and implementation of binding predictor
                                (version 1.0 2022/06/05)

1 Description    
    DisoBindPred is an accurate sequence-based predictor for identifying disordered binding residues and disordered protein-binding residues
	
2 Installation
    DisoBindPred is built on Python3.
    We highly recommend to use a virtual environment for the installation of DisoBindPred.

2.1 Create an virth environment
	Note: Before create the virth enviroment, make sure conda is installed. Otherwise, install Miniconda (https://conda.io/en/latest/miniconda.html) for conda.     
		e.g., On linux system, install Miniconda using 
		$ sh Miniconda3-py38_4.12.0-Linux-x86_64.sh -b
	        $ ~/miniconda3/bin/conda init	

    A virtual environment can be created and (de)activated as follows by using conda(https://conda.io/docs/):
	# create the virth environment disorderbind
	$ conda create -n DisoBindPred  python=3.8.12
	# activate the virth environment DisoBindPred 
	$ source activate DisoBindPred 
		# install package for machine learning algorithms
		$ conda install scikit-learn==0.24.2 kares==2.4.0 numpy==1.18.5 tensorflow==2.3.0
	#deactivate after finish the installation of DisoBindPred 
	$source deactivate

2.2 Install tools and database

	Note: $path is the absolute path for DisoBindPred installation

	2.2.1 All required tools
	(1) Check tcsh-6.21.00 in $path/programs/tcsh-6.21.00 for usage of PsiPred. If need, please recompile in the given directory. 
	(2) Check jre1.8.0_192 in $path/programs/jre1.8.0_192 for usage of VSL2b. If need, please recompile in the given directory. 
	(3) Check PsiPred in $path/programs/PSIPRED. If need, please recompile in the given directory. 
	(4) Check ASAquick in $path/programs/ASAquick. If need, please recompile in the given directory. 
	(5) Check VSL2B in $path/programs/VSL2B. If need, please recompile in the given directory. 
	(6) Check mmseqs in $path/programs/mmseqs. If need, please recompile in the given directory.
	(7) Check hh-suite in $path/programs/hhblits. If need, please recompile in the given directory. 
	(8) Check ncbi-blast-2.13.0+ in $path/programs/ncbi-blast-2.13.0+. If need, please recompile in the given directory.
	
	2.2.2 Databases

	Note: Please prepare two databases that are used for mmseqs and hhblits in $path/database.  
		

	(1) Please download database Uniref90 from (ftp://ftp.ncbi.nlm.nih.gov/blast/db/).
		# save in the directory $path/database/mmseqs_db
		# create dataset for mmseqs 
			$ cd $path/database/mmseqs_db
			$ $path/programs/mmseqs/bin/mmseqs createdb uniref90.fasta uniref90

	(2) Please download database uniclust30_2018_08 from (http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz).
		# save in the directory $path/database
		$ tar xvfz uniclust30_2018_08_hhsuite.tar.gz
	

2.3 Compile the scripts for DisoBindPred 

	Note: $path is the absolute path for DisoBindPred installation
	
	(1) compile aln_ver2.sh
		$ cd $path
		$ chmod a+x aln_ver2.sh
	$chmod a+x coverage.sh
	
	(2) Compile cov231stats.c	
		$ cd $path
		$ gcc -O3 -ffast-math cov231stats.c -o cov231stats

	(3) change the directory of environment to the virth environment DisoBindPred 
	 	e.g. source ~/miniconda3/bin/activate ~/miniconda3/envs/DisoBindPred #change the environment directory


2.4 Run DisoBindPred 

	Note: $path is the absolute path for DisoBindPred installation. Input_file should be in FASTA format. 
		
	To run the predictor, please use the following command.
		$ python run_binding_predictor.py test_seq.txt $path test_output(where output will be saved)	
		Note: $path is where the DisoBindPred installed. 

2.5 The test input file (test_seq.txt) and outputs are saved in $path/test_output

