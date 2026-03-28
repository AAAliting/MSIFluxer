# -------------------------------------------------------------------------------
# Name:        mfapy example 3 13C-MFA of brest cancer (MCF-7) cells
#              Metabolic model and data used in this code is derived from Araki et al Mass Spectrometry 2018, 7, A0067.
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
# -------------------------------------------------------------------------------
import mfapy
import os, sys, time
import pandas as pd
import numpy as np

if __name__ == '__main__':

    results_df = pd.DataFrame(columns=["Method", "RSS1", "RSS2", "RSS3", "RSS4", "RSS5"])
    # Total 100 (=20*5) random initial flux vectors were prepared
    numofmodel = 20
    repeatn = 5
    #
    # Construction of metabolic model
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("MFA_model.txt",format="text")
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments, mode="normal")
    #
    # Configurations
    #
    model.set_configuration(callbacklevel=0)  
    model.set_configuration(iteration_max=10000)  
    model.set_configuration(ncpus=4)  # Number of local CPUs for joblib
    #
    # Batch setting of constrains from  file
    #
    state_dic = model.load_states("SC_status.csv", format='csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Set isotope labelling of carbon sources
    #
    carbon_source1 = model.generate_carbon_source_template()
    carbon_source1.set_all_isotopomers('SubsCO2', [0.99, 0.01])
    carbon_source1.set_each_isotopomer('SubsGlc', {'#000000': 0.47, '#111111': 0.53}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsGln', {'#00000': 1.}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsSer', {'#000': 1.}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsAla', {'#000': 1.}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsLac', {'#000': 0.47, '#111': 0.53}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsR5P', {'#00000': 1.}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsAcCOA', {'#00': 1.}, correction='yes')
    carbon_source1.set_each_isotopomer('SubsGly', {'#000': 1.}, correction='yes')
    #
    # Load measured MDV data
    #
    mdv_observed1 = model.load_mdv_data('SC_MDV.txt', )
    mdv_observed1.set_std(0.015, method='absolute')
    #
    # Set experiment
    #
    model.set_experiment('ex1', mdv_observed1, carbon_source1)
    #
    #
    print("Generation of initial states using parallel processing using joblib") 
    flux_initial_array = []
    rss_initial_array = []
    for j in range(repeatn):
        state, flux_initial = model.generate_initial_states(numofmodel * 10, numofmodel, method="parallel")
        flux_initial_array.append(flux_initial)
        rss_initial_array.extend(model.calc_rss(flux_initial))
    temparray = sorted(rss_initial_array)
    datan = len(temparray) - 1
    print("initial", datan,"RSS:{0:>8.2f}  RSS:{1:>8.2f} RSS:{2:>8.2f} RSS:{3:>8.2f} RSS:{4:>8.2f}".format(temparray[0],temparray[int(datan * 0.05)], temparray[int(datan * 0.5)],temparray[int(datan * 0.95)],temparray[datan]))

    model.set_configuration(callbacklevel=0)  
    #
    # Fitting by deep optimizer
    #
    RSS_dict = {}

    all_methods = ["SLSQP", "COBYLA", "LN_COBYLA", "LN_BOBYQA", "LN_PRAXIS", "LN_SBPLX", "LN_NELDERMEAD", "LN_NEWUOA", "GN_CRS2_LM", "GN_DIRECT_L", "GN_AGS", "GN_IRES", "GN_ESCH"]  #

    for method in all_methods:
        RSS_dict[method] = []
        for j in range(repeatn):
            start = time.time()
            flux_opt1 = flux_initial_array[j]
            for i in range(1000):
                state, RSS_bestfit, flux_opt1 = model.fitting_flux(method=method, flux=flux_opt1)
                if len(flux_opt1) == 0:
                    break
                pvalue, rss_thres = model.goodness_of_fit(flux_opt1[0], alpha=0.05)
                now = time.time()
                if now - start > 60:
                    break
            RSS_dict[method].extend(RSS_bestfit)
        RSS_dict[method].sort()
        temparray = RSS_dict[method]
        if len(temparray) == 0:
            continue
        datan = len(temparray) - 1

        print(method, datan,
              "RSS:{0:>8.2f}  RSS:{1:>8.2f} RSS:{2:>8.2f} RSS:{3:>8.2f} RSS:{4:>8.2f}".format(temparray[0], temparray[
                  int(datan * 0.05)], temparray[int(datan * 0.5)], temparray[int(datan * 0.95)], temparray[datan]))

        # DataFrame
        results_df = pd.concat([results_df, pd.DataFrame({"Method": [method],
                                                          "RSS1": [temparray[0]],
                                                          "RSS2": [temparray[int(datan * 0.05)]],
                                                          "RSS3": [temparray[int(datan * 0.5)]],
                                                          "RSS4": [temparray[int(datan * 0.95)]],
                                                          "RSS5": [temparray[datan]]})], ignore_index=True)

    # write result to Excel
    results_df.to_excel("method_comparison_results.xlsx", index=False)

# read
    df = pd.read_excel("method_comparison_results.xlsx")

    # name
    methods = df['Method']

    # get median
    median_results = df.iloc[:, 1:].median(axis=1)

 
    median_results_top8 = median_results.iloc[:8]
    median_results_last5 = median_results.iloc[8:]

    
    min_median_method_top8 = methods.iloc[median_results_top8.idxmin()]

   
    min_median_method_last5 = methods.iloc[median_results_last5.idxmin()]

    print(f"local ：{min_median_method_top8}")
    print(f"global：{min_median_method_last5}")
    
    with open("method_comparison_results.txt", "w") as f:
        if min_median_method_top8 is not None:
            f.write(f"min_median_method1: {min_median_method_top8}\n")
        else:
            f.write("none\n")

        if min_median_method_last5 is not None:
            f.write(f"min_median_method2: {min_median_method_last5}")
        else:
            f.write("none")
# MFA analysis
start = time.time()

if __name__ == '__main__':
    #
    # Construction of metabolic model
    start = time.time()
    print(start)
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("MFA_model.txt")
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel=4)
    model.set_configuration(iteration_max=1000000)  # Maximal iternations in optimization
    model.set_configuration(
        number_of_repeat=3)  # Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus=5)  # Number of local CPUs for Parallel python
    #
    # Loading metabolite state from text file to constrain flux vector
    #
    state_dic = model.load_states("Status_Folder/MFA_status.csv", format='csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    carbon_ratio = pd.read_csv('MFA_Carbon_Source.csv', index_col=0)
    # read method
    with open("method_comparison_results.txt", "r") as f:
        min_median_method1 = f.readline().split(":")[1].strip()
        min_median_method2 = f.readline().split(":")[1].strip()

    source_dir1 = 'MDV_MFA/'
    dest_dir1 = 'Results_MFA/'
    files1 = os.listdir(source_dir1)

    for file in files1:
        print(file)
        result_path = dest_dir1 + file
        C12 = carbon_ratio.loc['C12', file]
        C13 = carbon_ratio.loc['C13', file]
        carbon_source1 = model.generate_carbon_source_template()
        carbon_source1.set_all_isotopomers('SubsCO2', [0.99, 0.01])
        carbon_source1.set_each_isotopomer('SubsGlc', {'#000000': C12, '#111111': C13}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsGln', {'#00000': 1.}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsSer', {'#000': 1.}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsAla', {'#000': 1.}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsLac', {'#000': C12, '#111': C13}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsR5P', {'#00000': 1.}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsAcCOA', {'#00': 1.}, correction='yes')
        carbon_source1.set_each_isotopomer('SubsGly', {'#000': 1.}, correction='yes')
        mdv_observed1 = model.load_mdv_data(os.path.join(source_dir1 + file),format='txt')
        mdv_observed1.set_std(0.015, method='absolute')
        model.set_experiment('ex1', mdv_observed1, carbon_source1)  #
        method_global = min_median_method2
        # methods = ["SLSQP", "COBYLA", "LN_COBYLA", "LN_BOBYQA", "LN_NEWUOA", "LN_PRAXIS", "LN_SBPLX",
        #           "LN_NELDERMEAD", "GN_CRS2_LM", "deep"]
        method_local = min_median_method1
        state1, flux_opt1 = model.generate_initial_states(2000, 5, method="parallel")
        results = [('template', flux_opt1[0])]

        state, RSS_bestfit, flux_opt1 = model.fitting_flux(method=method_global, flux=flux_opt1)
        for i in range(5):
            results.append((method_global, flux_opt1[i]))

        state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method=method_local, flux=flux_opt1)
        for i in range(5):
            pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp[i], alpha=0.05)
            results.append((method_local, flux_opt_slsqp[i]))

        model.show_results(results, filename=result_path + '.csv', format="csv")

end = time.time()
duration = end - start
print(duration)
