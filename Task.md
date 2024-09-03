# Required Tasks

## Final goals
#### Input data:
1. Single output
2. Larger number of Single output from monitoring, not sync, not same length, not same excitation ($N \geq 3000$ )
3. Focus on signal is 1OG, Z dirs.
4. Menhir data (Vel)
5. No excitation (input data), and should not assume they are white noises.

#### Output result:
1. System identification with uncertainty. ideal: also extract natural frequency and dampin ratio.
2. FRFs with uncertainty with ROM, can select the number of DOF formed.
3. Investigate the potential input excitation for each event (might need clustering first)
4. Surrogate model that can predict the FRFs (but what is input?)

#### Framework
![Overall framework](/MDSI_FWt.svg)


## TODO Task List: 

1. Validate the proposed EMA-FRF technique with OMA-COV Acc mesuremnt results in scatter plot.
2. Further tasks

## Finish and Conclusion:

#### Visualization 
1. Visualization: [Main_Visualized_all_Menhir_data.m](/MATLAB/Main/Main_Visualized_all_Menhir_data.m)
   -  Visualized All fo the Menhir data for 1G Z dirs
  
#### EMA testing 
2. EMA testing: [Main_EMA_Acc_Measure.m](/MATLAB/Main/Main_EMA_Acc_Measure.m)
    - Data is from Acc measurement.
    - Find the FRFs of recorded point 1G, 15 records in Acc.
    - Acc data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Ouput signal: channel 9, Input signal: channel 19, Using PSD method.
    - Conclusion: Use this as the validation result in the future. Potential Natural frequency are [12.5, 16.5, 19, 62.5] (Hz).

3. EMA testing: [Main_EMA_Menhr_Measure_1OG_Z.m](/MATLAB/Main/Main_EMA_Menhr_Measure_1OG_Z.m)
    - Data is from Menhir database. Only use the data where sensor in 1G and 2G are activated simultaneously.
    - Find the FRFs of recorded point 1G, in Vel.
    - Data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Ouput signal: Z_mm_s, Input signal: channel 19, Using FFT method.
    - Conclusion: Potential Natural frequency are [11.96, 19] (Hz).
  
4. EMA testing: [Main_ComparisonManhrAcc_Measure_Hammer_FRF.m](/MATLAB/Main/Main_ComparisonManhrAcc_Measure_Hammer_FRF.m)
    - Compare the data of Menhir (Vel) and Acc measurement (ACC) one by one.
    - Data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Only focus on the Hz smaller than 50 Hz.
    - The FRFs are all noramlized for comparison.
    - Conclusion: After normalization, the FRFs from ACC and Menhir are similar, where the Potential Natural frequency are appr. [12, 19] (Hz),however, the damping ratio are unclear.

#### OMA testing 

5. OMA testing: [Main_ComparisonManhrAcc_Measure_1OG_Z.m](/MATLAB/Main/Main_ComparisonManhrAcc_Measure_1OG_Z.m)
    - Compare the Data from Menhir (Vel) and Acc measurment (ACC) using [SINGLE-SENSOR system identification method](/MATLAB/Lib/ECheynet-modalID_singleSensor-b67ee9f/).
    - Data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Focus on the data in the first floor, Z dirs.
    - Conclusion: From both Menhir data and Acc data, their natrual frequencies detected are similar, however the damping ratio are a little bit diverge, but overall, they shows the similar results.
    - Concern: eventhough the dampring raito deteced are in similar range, however this range might be incorrect.

6. OMA testing: [Main_ComparisonManhrAcc_Measure_2OG_Z.m](/MATLAB/Main/Main_ComparisonManhrAcc_Measure_2OG_Z.m)
    - Compare the Data from Menhir (Vel) and Acc measurment (ACC) using [SINGLE-SENSOR system identification method](/MATLAB/Lib/ECheynet-modalID_singleSensor-b67ee9f/).
    - Focus on the data in the second floor, Z dirs.
    - Data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Conclusion: From both Menhir data and Acc data, their natrual frequencies detected are similar, however the damping ratio are a little bit diverge, but overall, they shows the similar results.
    - Limitation: Here we only focus on the Freq range between [0, 30] Hz, if the detected range is larger than this, then the result looks different.
    - Concern: eventhough the dampring raito deteced are in similar range, however this range might be incorrect.
  
7. OMA testing: [Main_OMA_SSI_COV_Acc_Measure_1OG_Z.m](/MATLAB/Main/Main_OMA_SSI_COV_Acc_Measure_1OG_Z.m)
    - System identification using the [SSI-COV mathod](/MATLAB/Lib/ECheynet-modalID_singleSensor-b67ee9f/)
    - Focus on the data in the first floor, Z dirs.
    - Data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Focus on the ACC data, since we have serveral Acc data recorded simultaneously from the same slab ar 1 O.G Z dir. -> suitable for SSI-COV method.
    - Conclusion: Similar natrual frequency as using SINGLE-SENSOR system identification method, however different damping ratio are reported.
    - Limitation: We need serveral ouput data from the structure at the same time (MultiOutput problem).
    - Assumption: Using SSI-COV technique might give us better results. 
  
8. OMA testing: [Main_OMA_SSI_COV_Acc_Measure_2OG_Z.m](/MATLAB/Main/Main_OMA_SSI_COV_Acc_Measure_2OG_Z.m)
    - System identification using the [SSI-COV mathod](/MATLAB/Lib/ECheynet-modalID_singleSensor-b67ee9f/)
    - Focus on the data in the second floor, Z dirs.
    - Data is recorded by Hammer Measurement, experimental date: 12.03.2024.
    - Focus on the ACC data, since we have serveral Acc data recorded simultaneously from the same slab ar 2 O.G Z dir. -> suitable for SSI-COV method.
    - Conclusion: Similar natrual frequency as using SINGLE-SENSOR system identification method, however different damping ratio are reported.
    - Limitation: We need serveral ouput data from the structure at the same time (MultiOutput problem). 
    - Assumption: Using SSI-COV technique might give us better results.

#### Forming FRFs

9. Forming FRFs: [Main_EMA_OMA_FRF_Acc_Measure_compare.m](/MATLAB/Main/Main_EMA_OMA_FRF_Acc_Measure_compare.m)
    - Compare the FRF of Acc measuremnt and FRF formed by superposition of 2 DOF system. 
    - Using OMA technique.
    - Using [SINGLE-SENSOR system identification method](/MATLAB/Lib/ECheynet-modalID_singleSensor-b67ee9f/) technique to extract the first two natural frequencies and damping ratios, after that, forming the FRF of SDOF system based on each natural frequencies and damping ratio, then superposite them. 
    - Conclusion: similar natural frequencies but differnet amplitude and coverd energey, meaning the damping ratio detected are different (bad estimation of damping ratio)

10. Forming FRFs: [Main_EMA_SI_FRF_form_ACC_Measure_1OG_Z.m](/MATLAB/Main/Main_EMA_SI_FRF_form_ACC_Measure_1OG_Z.m)
    - Using EMA technique.
    - Using rational based polynomial. 
    - Extracted natrual frequencies and damping ratios from fitting rational based polynomial. 
    - Form ROM FRF, specifing the number of DOF . 
    - Self-written method, didn't validate yet. 
    - Conclusion: similar FRFs as EMA results, however the amplitude and damping ratio sometimes are still not fit well. 
  

#### Clustering (03.09.2024)
11. Clustering using hierarchical cluster(HC): [Main_clustering_test_HC_2.m](/MATLAB/Main/Main_clustering_test_HC_2.m)
    - Focus on the Menhir data in z dir, 1 O.G. $N$ = 1528
    - Using hierarchical cluster methd (linkage)
    - Better than DBSCAN in mu opinion. 
    - Ignore the cluster that the number of events include is smaller than 50 (${N_{remain} = 1365}$  ).
    - Based on move-mean normalized signal ([0,1]), meaning that here I only care about the shape, however the aplitudes are not considered well here. 
    - Default: 15 clusters for no reason -> might need to find a better way to distiguish the number of cluster we want. 



##  Questions and methods

#### For clustering: 
- Identify the singal by their enegry (PSD) and Pattern (PSD noramlized with energy)
- How to extract the pattern? Extact top three peaks and damping ratio? 
- How to denoise? using PSD? 
- What if my sinal have different size? Transfer to Frequency domain to eliminate the time info? 
- How to cluster? using K-mean? HC? DBSCAN? or other method? 
- What is the expected result? Clustering by Ampitude and pattern? What if same pattern but different amplitude? 
- System is stationary, however the input singal and output singal should be non-statinoary (need SFFT or Wavelet?)




