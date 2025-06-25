# precision_cancer_metabolism_targeting #

### How to install
---
It works on Linux operating system and has been tested on Ubuntu 16.04.

1. Clone the repository
```
git clone https://github.com/kaist-sbml/precision_cancer_metabolism_targeting.git
```
2. Create and activate a conda environment
```
conda env create -f environment.yaml
```

---
### Implementation
1. Create files for flux values, context-specific GEMs and reaction weights. 
Every input file must follow the required naming convention: <sample_name>_flux.csv, <sample_name>_GEM.xml, and <sample_name>_reaction_weight.csv.

(See the example in example_input directory)

2. Run knock-down simulation. The following arguments are required: -o/--output_dir, -i/--input_dir
	
	#### Example
	```
	python run_KDsimulator.py -i ./example_input/ -o ./example_output/
	```

---
### Reference
Publication will be added here.