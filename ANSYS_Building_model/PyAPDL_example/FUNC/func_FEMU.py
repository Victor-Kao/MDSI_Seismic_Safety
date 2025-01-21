from Simulation_PyAPDL import simulation_PyAPDL
from sklearn.cluster import DBSCAN
from scipy.stats import mode
from scipy.integrate import simps
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import json
import os
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import seaborn as sns
from scipy.stats import spearmanr
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol
from mpl_toolkits.mplot3d import Axes3D


def get_MAC(mode_shape_model,mode_shape_exp):
    mode_shape_model = mode_shape_model.reshape(1, -1)
    mode_shape_exp = mode_shape_exp.reshape(1, -1)
    return np.power(np.dot(mode_shape_exp,np.transpose(mode_shape_model)),2)[0][0]/(np.dot(mode_shape_model,np.transpose(mode_shape_model))[0][0]*np.dot(mode_shape_exp,np.transpose(mode_shape_exp)))[0][0]

def fit_gaussian_kernel(f,f_n,zeta):
    # Compute the standard deviation (sigma) based on the damping ratio
    band_width_ratio = 1.1
    sigma = band_width_ratio*(zeta * f_n) / np.sqrt(2 * np.log(2))
    return np.exp(-0.5 * ((f - f_n) / sigma)**2)

def find_peaks_SDOFsup(num_peak,natrual_freq_arr,frf_freq_arr,targer_FRF,zeta, plot_ = False):
    frf_freq_arr = frf_freq_arr.reshape(1,-1)[0]
    targer_FRF = targer_FRF.reshape(1,-1)[0]
   
    peak_list = np.zeros(num_peak)
    peak_freq_list = np.zeros(num_peak)
    peak_energy_list = np.zeros(num_peak)
    kernel = np.zeros(len(frf_freq_arr))
    kernel_current_best = np.zeros(len(frf_freq_arr))
    
    for i_peak in range(num_peak):
        energy_remain_arr = np.ones(len(natrual_freq_arr))
 
        for i_freq_n in range(len(natrual_freq_arr)):#
            freq_target = natrual_freq_arr[i_freq_n]
            ampl = np.interp(freq_target, frf_freq_arr, targer_FRF)
            kernel = np.maximum(kernel_current_best, ampl*fit_gaussian_kernel(frf_freq_arr,freq_target,zeta))

            signal_remain = targer_FRF - kernel
            energy_remain_arr[i_freq_n] = simps(signal_remain**2, frf_freq_arr)

        peak_energy_list[i_peak] = np.min(energy_remain_arr)
        peak_pos = np.argmin(energy_remain_arr)
        peak_list[i_peak] =  peak_pos
        peak_freq_list[i_peak]  = natrual_freq_arr[peak_pos]
        
        freq_current_best = natrual_freq_arr[peak_pos]
        ampl_current_best= np.interp(freq_current_best, frf_freq_arr, targer_FRF)
        kernel_current_best = np.maximum(kernel_current_best, ampl_current_best*fit_gaussian_kernel(frf_freq_arr,freq_current_best,zeta))

    if plot_:
        plt.plot(frf_freq_arr, targer_FRF,linestyle=":")
        plt.plot(frf_freq_arr, kernel_current_best,linestyle=":")
        plt.plot(frf_freq_arr, targer_FRF -kernel_current_best ,linestyle=":")
        plt.show()
    
    
    return np.sort(peak_list),np.sort(peak_freq_list)

def extract_mode_shape_vector(f_n_arr, frf_arr):
    ## Get mode shape from FRF
    i_model_1OG = np.array([ frf_arr['disp_ch9']['imag'].reshape(1,-1),
                            frf_arr['disp_ch10']['imag'].reshape(1,-1),
                            frf_arr['disp_ch11']['imag'].reshape(1,-1),
                            frf_arr['disp_ch12']['imag'].reshape(1,-1)])
    i_model_1OG_mat = np.vstack(i_model_1OG)
    max_i_1 = np.max(abs(i_model_1OG_mat))
    i_model_1OG_norm = i_model_1OG_mat/max_i_1

    i_model_2OG = np.array([ frf_arr['disp_ch3']['imag'].reshape(1,-1),
                            frf_arr['disp_ch13']['imag'].reshape(1,-1),
                            frf_arr['disp_ch14']['imag'].reshape(1,-1),
                            frf_arr['disp_ch15']['imag'].reshape(1,-1)])
    i_model_2OG_mat = np.vstack(i_model_2OG)
    max_i_2 = np.max(abs(i_model_2OG_mat))
    i_model_2OG_norm = i_model_2OG_mat/max_i_2

    f = frf_arr['disp_ch9']['freq'].reshape(1,-1)
    mode_shape_vector = np.zeros([len(f_n_arr),8])
    mode_freq_vector = np.zeros([len(f_n_arr)])

    i_mode = 0
    ampl_ratio = max_i_1/max_i_2
    for f_i_mode in f_n_arr:
        i_shape = 0
        mode_freq_vector[i_mode] = f_i_mode
        for i_mode_shape in range(4):
            mode_shape_1OG = np.interp(f_i_mode, f[0], i_model_1OG_norm[i_mode_shape,:])
            mode_shape_vector[i_mode,i_shape] = mode_shape_1OG
            i_shape = i_shape +1
            
        for i_mode_shape in range(4):
            mode_shape_2OG = np.interp(f_i_mode, f[0], i_model_2OG_norm[i_mode_shape,:]*ampl_ratio)
            mode_shape_vector[i_mode,i_shape] = mode_shape_2OG 
            i_shape = i_shape +1
        
        i_mode = i_mode +1

    return mode_freq_vector, mode_shape_vector

def mean_value_filted(data,std_dev_thres = 1):
    # Calculate the mean and standard deviation of the data
    mean = np.mean(data)
    std_dev = np.std(data)
    # Filter out outliers
    filtered_data = data[np.abs(data - mean) <= std_dev_thres * std_dev]

    if len(filtered_data) == 0:
        filtered_data = data
    # Calculate the mean of the remaining data
    return np.mean(filtered_data)

def find_max_cluster(data):
    # Reshape for clustering
    data_reshaped = data.reshape(-1, 1)
    # Apply DBSCAN for clustering
    db = DBSCAN(eps=1, min_samples=2).fit(data_reshaped)
    # Get cluster labels
    labels = db.labels_
    # Identify the largest cluster
    largest_cluster = mode(labels[labels != -1])[0][0]  # Exclude noise (-1)
    # Get values in the largest cluster
    return data[labels == largest_cluster]


## MMI :Improved finite element model updating of a full-scale steel bridge using sensitivity analysis　(Bjørn T. Svendsen, 2021)

def get_MMI(f_n_model,f_n_exp, mac, f_n_ratio=0.5):
    return(1-f_n_ratio)*mac - f_n_ratio* abs(f_n_exp-f_n_model)/f_n_exp

def split_data_limit_corr(input_X, input_Y, max_corr = 0.2, num_trial = 2000, test_size_ = 0.1, random_seed = None):
    for i_test in range(num_trial):
        bool_high_corr = False
        # Split into training and test sets
        if random_seed is None:
            seed = np.random.randint(0, 1000000)
            X_train, X_test, y_train, y_test = train_test_split(input_X, input_Y, test_size=test_size_, random_state=seed)

            # Compute Pearson and Spearman correlations for each feature
            correlations = []
            for i in range(input_X.shape[1]):
                pearson_corr = np.corrcoef(X_train[:, i][:X_test.shape[0]], X_test[:, i])[0, 1]
                spearman_corr, _ = spearmanr(X_train[:, i][:X_test.shape[0]], X_test[:, i])
                correlations.append((f"Feature {i+1}", pearson_corr, spearman_corr))
                if np.abs(pearson_corr) >= max_corr:
                    bool_high_corr = True
                    break

            if bool_high_corr:
                if i_test == num_trial-1:
                    print("Failed split data")
                continue
            else:
                print("Succesfully split data")
                print(f"Random seed {seed}")
                df = pd.DataFrame(correlations, columns=["Feature", "Pearson Corr", "Spearman Corr"])
                print(df)
                return X_train, X_test, y_train, y_test
        else:
            X_train, X_test, y_train, y_test = train_test_split(input_X, input_Y, test_size=test_size_, random_state = random_seed)
            return X_train, X_test, y_train, y_test
        
def sobol_GSA(gpr_model,num_var, sample_size = 2048, var_name = [], plot_ = False):
    # Define the problem (bounds of the input parameters)
    if var_name == []:
        problem = {
            'num_vars': 19,  # Number of parameters
            'names': ['param_{}'.format(i) for i in range(1, num_var+1)],  # Parameter names
            'bounds': [(0, 1)] * 19  # Parameter bounds, here [0, 1] for each
        }
    else:
        problem = {
            'num_vars': num_var,  # Number of parameters
            'names': [var_name[i] for i in range(0, num_var)],  # Parameter names
            'bounds': [(0, 1)] * num_var  # Parameter bounds, here [0, 1] for each
        }
    # Generate Sobol samples
    param_values = sobol_sample.sample(problem, sample_size)  
    Y = gpr_model.predict(param_values)
    # Perform Sobol analysis on the GP model's predictions
    Si = sobol.analyze(problem, Y, print_to_console=True)

    if plot_:
        # Create a figure and axis object
        fig, ax1 = plt.subplots(figsize=(12, 6))

        # Plot First-order Sobol indices (S1) on the left y-axis
        width = 0.35  # Width of the bars
        ax1.bar(var_name, Si['S1'], color='blue', label='First-order Sobol', alpha=0.7, width=width)
        ax1.set_xlabel('Parameters')
        ax1.set_ylabel('First-order Sobol indices', color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')
        ax1.set_xticklabels(var_name, rotation=90)

        # Create the second y-axis for Total Sobol indices (ST)
        ax2 = ax1.twinx()
        ax2.bar(var_name, Si['ST'], color='lightcoral', label='Total Sobol', alpha=0.7, width=width, align='edge')
        ax2.set_ylabel('Total Sobol indices', color='lightcoral')
        ax2.tick_params(axis='y', labelcolor='lightcoral')

        # Align the bars by adjusting the y-limits to the same range
        y_min = min(min(Si['S1']), min(Si['ST']))
        y_max = max(max(Si['S1']), max(Si['ST']))
        ax1.set_ylim([y_min-0.05, y_max+0.1])
        ax2.set_ylim([y_min-0.05, y_max+0.1])
        # Set title and adjust layout
        plt.title('First-order and Total Sobol indices (Sensitivity Analysis)')
        plt.tight_layout()

        # Show the plot
        plt.show()
        
    return Si


def visualizeSurrogate3D(gpr_model, num_var, var1, var2,input_name = []):
    # Define the indices of the two variables you want to vary
    var1_index = var1  # Index for the 5th variable (Python is zero-indexed)
    var2_index = var2  # Index for the 10th variable

    if input_name == []:
        input_name = []
        for i in range(num_var):
            input_name.append(i)
    else:
        input_name


    # Calculate the mean of the remaining 17 variables
    mean_vars = 0.5*np.ones(num_var)


    # Generate a grid of values for var1 and var2 (these will vary, others are fixed)
    var1_grid, var2_grid = np.meshgrid(np.linspace(0, 1, 50),
                                    np.linspace(0, 1, 50))

    # Flatten the grids for var1 and var2 to create 1D arrays
    var1_flat = var1_grid.ravel()
    var2_flat = var2_grid.ravel()

    # Create a new input matrix for predictions
    # Stack the selected variables and tile the other variables' mean values
    X_grid = np.column_stack([  np.tile(mean_vars, (len(var1_flat), 1))])
    X_grid[:,var1_index] = var1_flat
    X_grid[:,var2_index] = var2_flat
    # Generate predictions for the grid (assuming gpr is your trained model)
    y_pred_grid = gpr_model.predict(X_grid).reshape(var1_grid.shape)

    # Now, y_pred_grid has shape (50, 50), matching the grid's shape

    # Create a 3D figure
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot

    surf = ax.plot_surface(var1_grid, var2_grid, y_pred_grid, cmap='viridis', alpha=0.8)
    #ax.scatter(var1, var2, y_test, c='red', label='Test Results', s=50)

    # Add labels and title
    ax.set_xlabel(f'X{var1_index + 1}: {input_name[var1_index]}')
    ax.set_ylabel(f'X{var2_index + 1}: {input_name[var2_index]}')
    ax.set_zlabel('Predicted Y')

    ax.set_title('3D Surface Plot of surrogate model')
    # Add color bar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)

    # Show the plot
    plt.show()


def generate_mask(frequencies, ranges_and_factors):
    """
    Generate a mask for a signal based on specified frequency ranges and factors.
    
    Parameters:
        frequencies (numpy.ndarray): The frequency array.
        ranges_and_factors (list of tuples): A list of tuples, where each tuple contains:
            - A tuple specifying the frequency range (start, end).
            - A multiplicative factor for the range.
            
    Returns:
        numpy.ndarray: The generated mask.
    """
    mask = np.ones_like(frequencies)  # Initialize mask with 1
    for freq_range, factor in ranges_and_factors:
        start, end = freq_range
        mask[(frequencies >= start) & (frequencies <= end)] = factor
    return mask