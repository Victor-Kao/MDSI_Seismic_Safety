import csv
import numpy as np

def ExtractDictfromMenhirCsv(Menhir_path):

    # Initialize the dictionary with keys and empty lists
    data_dict = {
        'time': [],
        'X': [],
        'Y': [],
        'Z': []
    }

    # Open and read the CSV file
    with open(Menhir_path, mode='r') as file:
        csv_reader = csv.reader(file, delimiter=';')
        
        # Skip the first row (assumed to be a header or non-data row)
        next(csv_reader, None)
        
        time_start = 0.00
        time_int  = 0.001
        time_step = time_start

        i = 0
        # Iterate over the rows
        for row in csv_reader:
            time_step =  time_start + i*time_int
            #date_time = datetime.strptime(row[0], '%d.%m.%Y %H:%M:%S.%f')
            # Append data to the dictionary
            data_dict['time'].append(round(time_step,3))
            # Clean and convert numerical data
            data_dict['X'].append(float(row[1].replace(',', '.')))
            data_dict['Y'].append(float(row[2].replace(',', '.')))
            data_dict['Z'].append(float(row[3].replace(',', '.')))
            i = i +1

        # Convert lists to NumPy arrays
        for key in data_dict:
            data_dict[key] = np.array(data_dict[key])
        
        return data_dict