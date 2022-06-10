#!/bin/bash

#activate the environment for binding predictor (please change to where the virth environment inplemented)
source ~/miniconda3/bin/activate ~/miniconda3/envs/disbindpredict

python run_prediction.py $1 $2
python form_input_set.py $1 $2
python load_model.py $1 $2
python generate_output.py $1 $2