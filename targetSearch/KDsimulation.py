from cobra.io import read_sbml_model
import copy
import numpy as np
import pandas as pd
import math
import time
import warnings
import logging
import os
import glob
from targetSearch import LAD
from targetSearch import Simulator
import multiprocessing


class KDSimulator:
    def __init__(self, context_specific_cobra_model_file, reaction_weight_file, flux_file, target_mets):
        self.reaction_weights = pd.read_csv(reaction_weight_file, index_col=0).iloc[:,0].to_dict()
        self.context_specific_cobra_model = read_sbml_model(context_specific_cobra_model_file)
        self.flux_df = pd.read_csv(flux_file, index_col=0, header=0, names=["default"])
        sample_id = os.path.basename(flux_file)[:-4]
        stoichiometry_dict = {}
        for rxn in self.context_specific_cobra_model.reactions:
            stoichiometry_dict[rxn.id] = {}
            for met in rxn.metabolites:
                stoichiometry_dict[rxn.id][met.id] = abs(rxn.metabolites[met])

        self.stoichiometry_df = pd.DataFrame.from_dict(stoichiometry_dict).fillna(0)
        self.met_ids = [i.id for i in self.context_specific_cobra_model.metabolites]
#         target_mets = ['MNXM92_m', 'MNXM325_r', 'MNXM133_m', 'MNXM31869_e', 'MNXM12_c', 'MNXM105630_c', 'MNXM10_r', 'MNXM183_m', 'MNXM205_c']
        self.indices_of_traget_mets= [list(self.stoichiometry_df.index).index(x) for x in target_mets]
#         self.default_fluxsum = self.calculate_fluxsum(self.flux_df)['default']
#         pd.DataFrame.from_dict(self.default_fluxsum, orient='index', columns=['default']).to_csv("./fluxsum_TCGA_BLCA_All_210916/%s.csv"%sample_id)
        
    def calculate_fluxsum(self, flux_df):
        rxn_list = list(flux_df.index.values)
        error_rxn = []
        fluxsum_all = {}
        for each_col in flux_df:
            flux_val = flux_df.loc[self.stoichiometry_df.columns, each_col].fillna(0.0).abs()
            fluxsum_val = np.matmul(np.array(self.stoichiometry_df), np.array(flux_val), dtype=float)
            fluxsum_val = abs(fluxsum_val) * 0.5
        
            # for return the flux-sum of only target metabolites (metabolites consisting risk score formula) in context-specific GEM
#             fluxsum_all[each_col]={list(self.stoichiometry_df.index)[x]:fluxsum_val[x] for x in self.indices_of_traget_mets}
            
            # for return the flux-sum of all the metabolites in context-specific GEM
            fluxsum_all[each_col]={self.stoichiometry_df.index[x]:fluxsum_val[x] for x in range(self.stoichiometry_df.shape[0])} 
        
        return fluxsum_all
    
    def calculate_initial_flux(self, output_dir):
        obj = LAD.LAD()
        obj.load_cobra_model(self.context_specific_cobra_model)
        flux_constraints = {}
        flux_constraints['biomass_reaction'] = [0.01, 1000.0]   
        reaction_weight_info = {}
        model_reactions = [each_reaction.id for each_reaction in self.context_specific_cobra_model.reactions]
        for each_reaction in self.reaction_weights:
            if each_reaction in model_reactions:
                reaction_weight_info[each_reaction] = self.reaction_weights[each_reaction]/10000.0
                
        solution_status, objective_value, predicted_flux = obj.run_LP_fitting(opt_flux=reaction_weight_info, flux_constraints=flux_constraints)
        flux_df = pd.DataFrame.from_dict(predicted_flux, orient='index', columns=['flux'])
        flux_df.to_csv(output_dir)
        print("Saved as %s"%output_dir)
        
    def calculate_KDfluxsum(self, output_dir, bins=10, minimum_biomass=0.01):
        start = time.time()
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
        obj = LAD.LAD()
        obj.load_cobra_model(self.context_specific_cobra_model)

        model_reactions = [each_reaction.id for each_reaction in self.context_specific_cobra_model.reactions]
        model_flux_carrying_reactions = [each_reaction for each_reaction in model_reactions if abs(self.flux_df.loc[each_reaction].values) > 1e-12]
        
#         model_GPR_reactions = [each_reaction.id for each_reaction in self.context_specific_cobra_model.reactions
#                                if each_reaction.gene_reaction_rule != '']
#         print(len(model_GPR_reactions))

        print(len(model_flux_carrying_reactions))
        predicted_KD_fluxsum_all = {}
#         for each_GPR_reaction in model_GPR_reactions:
        for each_GPR_reaction in model_flux_carrying_reactions:
#             if each_GPR_reaction in self.reaction_weights and self.reaction_weights[each_GPR_reaction]!=0:
                predicted_KD_flux_all = {}
                for i in range(1, bins+1): 
                    new_reaction_weight_info = {}
                    for each_reaction in self.reaction_weights:
                        if each_reaction in model_reactions:
                            new_reaction_weight_info[each_reaction] = self.reaction_weights[each_reaction]/10000.0
                            
                    flux_constraints = {}
                    flux_constraints['biomass_reaction'] = [minimum_biomass, 1000.0]                               
                    constrained_target_flux = self.flux_df.loc[each_GPR_reaction].values*(1-1/bins*i)
                    if constrained_target_flux > 0:
                        flux_constraints[each_GPR_reaction] = [0, constrained_target_flux]
                    elif constrained_target_flux < 0:
                        flux_constraints[each_GPR_reaction] = [constrained_target_flux, 0]
#                     new_reaction_weight_info[each_GPR_reaction] = new_reaction_weight_info[each_GPR_reaction]*(1-1/bins*i)
#                     print(self.reaction_weights[each_GPR_reaction], new_reaction_weight_info[each_GPR_reaction])
                    else:
                        flux_constraints[each_GPR_reaction] = [0, 0]

                    solution_status, objective_value, predicted_flux = obj.run_LP_fitting(opt_flux=new_reaction_weight_info, flux_constraints=flux_constraints)
                    if solution_status ==2:
                        predicted_KD_flux_all[str(i*100/bins)+'%'] = predicted_flux
#                     break ## added for only save 0%
                    else:
                        print(solution_status, "Infeasible solution for {0}//{1}//{2} KD".format(os.path.basename(output_dir), each_GPR_reaction, str(i*100/bins)+'%'))
                        

                flux_df = pd.DataFrame.from_dict(predicted_KD_flux_all)
                ### added for save KD flux on 20220922
                flux_df.to_csv(output_dir+'/%s_KD_flux.csv'%each_GPR_reaction)
                ###
#                 print(flux_df)
#                 flux_df.columns = [str(x*100/bins)+'%' for x in range(1, bins+1)] 
#                 flux_df.columns = [str(x*100/bins)+'%' for x in range(0, 1)] ## added for only save 0%
### From here
                if len(flux_df.columns)>0:
                    predicted_KD_fluxsum_all[each_GPR_reaction] = self.calculate_fluxsum(flux_df)

                    # save KD fluxsum for KD each reaction
                    if not os.path.isdir(output_dir):
                        os.mkdir(output_dir)
                    pd.DataFrame.from_dict(predicted_KD_fluxsum_all[each_GPR_reaction]).to_csv(output_dir+'/%s_KD_fluxsum.csv'%each_GPR_reaction)
                else:
                    continue
                
                if len(predicted_KD_fluxsum_all)%500==0:
                    print("%s reactions done"%len(predicted_KD_fluxsum_all))
                    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
### to here commenting on 20220922 
#                     break
                
#             else:
#                 continue
        logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))        
        return predicted_KD_fluxsum_all