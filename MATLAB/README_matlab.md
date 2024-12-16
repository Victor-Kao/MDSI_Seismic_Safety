## MATLAB works for EMA/ OMA/ System identification
This folder includes the matlab files for system identificaiton.
Here are five sub-folders, which are: 
1. [Main](/Main/): All the analysis and some of the required ```.mat``` files are in this folder
2. [Func](/Func/): Includes all the required self-written function
3. [Lib](/Lib/): The library downloaded online or from the MATLAB community.
4. [Surrogate_main](/Surrogate_main/): Includes files and data used for builing surrogate model, such as DOE random generator
5. [MATALB_FIG](/MATALB_FIG/): Includes the figure that stored from analysis, used for presentation.

In this README, only [Main](/Main/), [Func](/Func/) and [Lib](/Lib/) are intrduced. the other folders are not highly relevent to system identification or analysis.

Author: Wei-Teng Kao (ge32gak@mytum.de)

## [Main folder](/Main/) 
In this section, only the files with the name start with "Main_" will explained, other files such as "TESTING_<file_name>" and the ```mat``` file are ignored.

In this folder, ```ACC``` in file's name denotes the acceleration data form hammer testng in 12.03.2024, which the accelerometer report the singal is accleration, this data can be used for EMA and OMA. 

On the other hand, ```Manhr``` in the file's name represent that the data used are extracted by the Manhir Geophone, which has been installed and monitored the operating vibration for over 2 years. This data can bed used for OMA only as we don't have the information of the input excitation.

### EMA
1. [Main_EMA_ACC_all_MeasPoint_per_test.m](/Main/Main_EMA_ACC_all_MeasPoint_per_test.m): Using the function [Func_PSD_FRF_COH.m](/Func/Func_PSD_FRF_COH.m) to form the FRFs at 8 measurement points on the slab for $z$ dir. Each floor has 4 measurment points and there are two floors. By selecting the ID of testing, it can formed the FRF based on ACC data. Here we have 13 valided testing case.

2. [Main_EMA_Acc_Measure_tfestimate_func.m](/Main/Main_EMA_Acc_Measure_tfestimate_func.m): Using the in-built function, tfestimate to form the FRFs at 8 measurement points on the slab for $z$ dir. Each floor has 4 measurment points and there are two floors. By selecting the ID of testing, it can formed the FRF based on ACC data. Here we have 13 valided testing case.
   
3. [Main_EMA_ComparisonManhrAcc_Measure_1OG_Z.m](/Main/Main_EMA_ComparisonManhrAcc_Measure_1OG_Z.m): Comparison of the FRF of Accelerometer and the PSD of the menhir data for the hammer testing. During the hammer testing, both of the Accelerometer and Manfir device recieved the excitation signal. This file testing the similarity of FRF extracted by ACC data and output of Menhir data. 
   
4. [Main_EMA_ComparisonManhrAcc_Measure_Hammer_FRF.m](/Main/Main_EMA_ComparisonManhrAcc_Measure_Hammer_FRF.m): Comparison of the FRF of the hammer testing. During the hammer testing, both of the Accelerometer and Manfir device recieved the excitation signal. This file testing the similarity of FRFs extracted from both devices.
   
5.  [Main_EMA_Menhr_EG_excitation.m](/Main/Main_EMA_Menhr_EG_excitation.m):  Form the FRF from operating data, not form the Hammer testing data. Here we ony focus on the Menhir data that both the Menhir sensors in EG, 1OG and 2OG are activated at the same time. In this analysis, I assumed the excitation is the singal recived from the sensor at the EG and the output is the signal recived from the 1OG or 2OG.

### OMA
1. [Main_OMA_ComparisonManhrAcc_Measure_[$n^{th}$ floor]_[$X,Y,Z$ dirs]_SingleSensorIden.m](/Main/Main_OMA_ComparisonManhrAcc_Measure_1OG_Z_SingleSensorIden.m): Compare the ACC and the Manhir data based on OMA Single_Sensor_Identification(SSI) method. Since this method is quite sensitvie, therefore the focusing freqeuncy range have to be specified based on each case.
2. [Main_OMA_Menhr_Measure_[$n^{th}$ floor]_SingleSensorIden.m](/Main/Main_OMA_Menhr_Measure_1OG_SingleSensorIden.m): Using Single_Sensor_Identification(SSI) package to identify the system. There are the cases for 1OG, 2OG and all three directions in individual files with similuar name. 
3. [Main_OMA_SSI_COV_ACC_Measure_1OG_Z_all_testing.m](/Main/Main_OMA_SSI_COV_ACC_Measure_1OG_Z_all_testing.m): Peforming the OMA using SSI-COV for ACC measurement data from 13 valid events in Hammer testing. The algorithm input the sensor data from 8 channels and output the potential natural frequencies, damping ratio and mode shapes. The excitations are at 2OG according to Hammer testing plan.
4. [Main_OMA_SSI_COV_Menhir_Measure_Z_all_floor.m](/Main/Main_OMA_SSI_COV_Menhir_Measure_Z_all_floor.m): Peforming the OMA using SSI-COV for Manhir Velocity data from several operating events. The algorithm input the Manhir sensor data at 1OG and 2OG floor (only 2 totally), and output the potential natural frequencies, damping ratio and mode shapes. The excitations are assumed to be the operation vibration propagating through the ground.

### Clustering
1. [Main_clustering_test_DBSCAN.m](/Main/Main_clustering_test_DBSCAN.m): Testing the ability of DBSCAN clustering.
2. [Main_clustering_test_HC_1.m](/Main/Main_clustering_test_HC_1.m): Testing the ability of HC clustering, which in my case is better than DBSCAN
, however DBSCAN has more possibilty. 



## [Func folder](/Func/)
Only the system identification realted function will be introduced.

#### Analysis used:
1. [Func_PSD_FRF_COH.m](/Func/Func_PSD_FRF_COH.m): Using PSD and CPSD to compute the FRF when input the excitation and the output signal. The input signal should be acceleration in time series, but the output is in displacement.
2. [Func_PSD_FRF_COH_freq.m](/Func/Func_PSD_FRF_COH_freq.m): Compute the FRF when input the excitation and the output signal. The input signal and output signal should be in frequency domain. There is no tranformation to displacement or vibration in this function.
3. [Func_HC_clustering.m](/Func/Func_HC_clustering.m): Using for clustering the signals, using in-built Hierarchical Clustering method
4. [Func_FFT_half.m](/Func/Func_FFT_half.m): FFT but only output the positive-half of the output result.
5. [Func_FilterDesign.m](/Func/Func_FilterDesign.m): Generating the band-pass filter.
6. [Func_FilterDesign_highlow.m](/Func/Func_FilterDesign_highlow.m): Generating the high-pass/ low-pass filter.
7. [Func_Resample.m](/Func/Func_Resample.m): Resample the data, this is used for adjusting the sampling rate for further analysis.
8. [Func_ConfiPlot.m](/Func/Func_ConfiPlot.m): Plot the result with the confident range, the input can be any kind of data, here I used it for input FRFs from each individual events/testings.
#### File managing related:
Those are the func used for importing/ transfering the data to MATALB readable/ compatible format. The raw data are stored in specific database.
1. [Func_FindCsvFiles.m](/Func/Func_FindCsvFiles.m)
2. [Func_FindDateTime.m](/Func/Func_FindDateTime.m)
3. [Func_FindMatFiles.m](/Func/Func_FindMatFiles.m)
4. [Func_ImportMenhirData2Tab.m](/Func/Func_ImportMenhirData2Tab.m)


## [Lib folder](/Lib/)
1. [Single_Sensor_Identification(SSI)](/Lib/ECheynet-modalID_singleSensor-b67ee9f/): Operational modal analysis with a single accelerometer. Citation: https://www.mathworks.com/matlabcentral/fileexchange/86718-operational-modal-analysis-with-a-single-accelerometer
2. [SSI_COV](/Lib/ECheynet-SSICOV-82ce27a/): Operational modal analysis with automated SSI-COV algorithm. Citation: https://www.mathworks.com/matlabcentral/fileexchange/69030-operational-modal-analysis-with-automated-ssi-cov-algorithm
3. [PickPeaks](/Lib/pickpeaks/): Using for finding the peaks from the signal. Citation: https://www.mathworks.com/matlabcentral/fileexchange/42927-pickpeaks-v-select-display

