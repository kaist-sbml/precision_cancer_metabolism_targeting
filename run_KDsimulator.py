from targetSearch.KDsimulation import KDSimulator
import pandas as pd
import glob, os
import numpy as np
import multiprocessing

def run_KDSimulator(model, weight, flux, target_mets, output_dir):
    for i in range(len(model)):
        simul = KDSimulator(model[i], weight[i], flux[i], target_mets)
        simul.calculate_KDfluxsum(output_dir[i], bins=5, minimum_biomass=0.01)
#         simul.calculate_initial_flux(output_dir[i])

    
if __name__=="__main__":
    num_cores = 30
    target_mets = ['MNXM92_m', 'MNXM325_r', 'MNXM133_m', 'MNXM31869_e', 'MNXM12_c', 'MNXM105630_c', 'MNXM10_r', 'MNXM183_m', 'MNXM205_c']
    
    coxData = pd.read_csv("./for_analysis/final_TCGA_BLCA_dataset_for_Cox_regression_211022.csv", index_col=0)
    sampleLst = [x+'-01' for x in coxData.index if coxData.loc[x,'risk_groups']=='High risk']
    outputDir = "./KD_fluxsum_TCGA_BLCA_High_risk_220924"
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    
    modelPaths = []
    fluxPaths = []
    weightPaths = []
    outputPaths = []
    for p in glob.glob("../input/TcgaBLCARNASeq_*"):
        check = [m for m in glob.glob(p+"/condition1/reconstructed_models/*") if os.path.basename(m) in sampleLst]
#         check = [m for m in glob.glob(p+"/condition1/reconstructed_models/*")]
        if len(check) > 0:
            modelPaths += [m+"/functional_%s.xml"%os.path.basename(m) for m in check]
            fluxPaths += ["./initial_flux_210916/%s.csv"%os.path.basename(x) for x in check]
            weightPaths += [p+"/flux1/Reaction_weight_%s.csv"%os.path.basename(x) for x in check]
            outputPaths += [outputDir+"/%s"%os.path.basename(x) for x in check]

    split_model_files = np.array_split(modelPaths, num_cores)
    split_flux_files = np.array_split(fluxPaths, num_cores)
    split_rxn_weights = np.array_split(weightPaths, num_cores)
    split_output_dirs = np.array_split(outputPaths, num_cores)
    
    pool = multiprocessing.Pool(num_cores)
    # use starmap to send multiple arguments
    pool.starmap(run_KDSimulator, [(split_model_files[i], split_rxn_weights[i], split_flux_files[i], target_mets, split_output_dirs[i]) for i in range(num_cores)])
    #[(result, local_dict, model_paths) for model_paths in all_paths]
    pool.close() 
    pool.join()