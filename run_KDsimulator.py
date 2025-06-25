from targetSearch.KDsimulation import KDSimulator
import pandas as pd
import glob, os
import numpy as np
import argparse

def argument_parser():
    parser = argparse.ArgumentParser()    
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory')
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")

    return parser

def run_KDSimulator(model, weight, flux, target_mets, output_dir):
    for i in range(len(model)):
        simul = KDSimulator(model[i], weight[i], flux[i], target_mets)
        simul.calculate_KDfluxsum(output_dir[i], strength=80, minimum_biomass=0.01)

    
if __name__=="__main__":
    target_mets = ['MNXM92_m', 'MNXM325_r', 'MNXM133_m', 'MNXM31869_e', 'MNXM12_c', 'MNXM105630_c', 'MNXM10_r', 'MNXM183_m', 'MNXM205_c']
    
    parser = argument_parser()
    options = parser.parse_args()
    
    input_dir = options.input_dir
    output_dir = options.output_dir
    
    for sample_flux in glob.glob(input_dir+'/*_flux.csv'):
        sample_model = sample_flux.replace('flux.csv', 'GEM.xml')
        sample_rw = sample_flux.replace('flux', 'reaction_weight')
        output_dir_sample = output_dir+'/'+sample_flux.split('/')[-1].split('_flux')[0]
        if not os.path.isdir(output_dir_sample):
            os.mkdir(output_dir_sample)

    
        _simulator = KDSimulator(sample_model, sample_rw, sample_flux, target_mets)
        _ = _simulator.calculate_KDfluxsum(output_dir_sample,strength=80, minimum_biomass=0.01)
