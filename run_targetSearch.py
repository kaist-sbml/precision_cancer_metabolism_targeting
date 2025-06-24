from targetSearch.KDsimulation import KDSimulator
from targetSearch.evaluation import *
import glob, os
import pandas as pd
import multiprocessing

def argument_parser():
    parser = argparse.ArgumentParser()    
    parser.add_argument('-n', '--num_cores', default= 13, required=True, help="The number of using cores")
    parser.add_argument('-f', '--flux', required=True, help="Input flux file directory")
    parser.add_argument('-w', '--weight', required=True, help="Input reaction weight file directory")
    return parser

def scaling(minmax_info, fluxsum):
    for met in fluxsum:
#         print(met, fluxsum[met])
        fluxsum[met] = (fluxsum[met]-minmax_info.loc['min'][met])/(minmax_info.loc['max'][met]-minmax_info.loc['min'][met])
        fluxsum[met]*=10
#         print(met, fluxsum[met])
        
    return fluxsum
    
    

def run_KDsimulation(model, weight, flux):
    
    target_mets = ['MNXM92_m', 'MNXM325_r', 'MNXM133_m', 'MNXM31869_e', 'MNXM12_c', 'MNXM105630_c', 'MNXM10_r', 'MNXM183_m', 'MNXM205_c']
    minmax_info =pd.read_csv("./BLCA_target_metabolites_fluxsum_minmax_info.csv", index_col=0)
    test = KDSimulator(model, weight, flux, target_mets)

    coefs = pd.read_csv("./target_met_coefs.csv", index_col=0).iloc[:,0].to_dict()
    result = test.calculate_KDfluxsum(bins=10, minimum_biomass=0.01)
    default_fluxsum = scaling(minmax_info, test.default_fluxsum)
    Topology=getLinearLinkage(test.context_specific_cobra_model, target_mets)
    score_dict={}
    for target_rxn in result:
        score_dict[target_rxn]=[[],[]]
        for KD_percent in result[target_rxn]:
            print(KD_percent)
            KD_fluxsum = scaling(minmax_info, result[target_rxn][KD_percent])
            AdjScore, Score = evaluate(default_fluxsum, KD_fluxsum, target_rxn, coefs, Topology, target_mets)

            score_dict[target_rxn][0].append(AdjScore)
            score_dict[target_rxn][1].append(Score)

    return score_dict
    



if __name__=="__main__":
    parser = argument_parser()
    options = parser.parse_args()
    num_cores = options.num_cores
    split_flux_files = np.array_split(glob.glob(options.flux), num_cores)
    split_rxn_weights = [[re.sub("Flux_prediction","Reaction_weight", x) for x in lst] for lst in split_flux_files]
    fluxSumData = pd.read_csv("./BLCA/testfluxsumresult_BLCA.csv", index_col=0)
    
    pool = multiprocessing.Pool(num_cores)
    # use starmap to send multiple arguments
    pool.starmap(run_KDsimulation, split_flux_files, split_rxn_weights)
    #[(result, local_dict, model_paths) for model_paths in all_paths]
    pool.close() 

    pool.join()