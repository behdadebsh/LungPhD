import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from openpyxl import load_workbook


def read_points_from_excel(file_path, sheet_name):
    df = pd.read_excel(file_path, sheet_name=sheet_name, engine='openpyxl')
    t_coords = df.iloc[:, 0].values
    y_coords = df.iloc[:, 1].values
    return t_coords, y_coords


def compute_fft_coefficients(y_coords, num_coefficients):
    N = len(y_coords)

    # Compute the FFT
    fft_coeffs = np.fft.fft(y_coords)


    # Extract the A_0 coefficient separately (DC component)
    A0 = fft_coeffs[0].real / N

    # Initialize arrays to store A_n and B_n coefficients
    An = np.zeros(num_coefficients)
    Bn = np.zeros(num_coefficients)

    for n in range(1, num_coefficients + 1):
        An[n - 1] = 2 * fft_coeffs[n].real / N
        Bn[n - 1] = -2 * fft_coeffs[n].imag / N

    return A0, An, Bn


def plot_fourier_series(t_coords, y_coords, A0, An, Bn, num_coefficients):
    T = t_coords[-1] - t_coords[0]
    t_new = np.linspace(t_coords[0], t_coords[-1], 100)
    y_fourier = np.full_like(t_new, A0)  # Initialize with A0
    y_new = np.full_like(t_new, A0)  # Initialize with A0
    # num_coefficients = len(An)
    Dn = np.zeros(num_coefficients)
    phi = np.zeros(num_coefficients)

    for n in range(1, num_coefficients + 1):
        Dn[n - 1] = np.sqrt(An[n - 1] ** 2 + Bn[n - 1] ** 2)
        phi[n - 1] = np.arctan2(Bn[n - 1], An[n - 1])

    for n in range(1, num_coefficients + 1):
        y_fourier += An[n - 1] * np.cos(2 * np.pi * n  * t_new / T) + Bn[n - 1] * np.sin(2 * np.pi * n * t_new / T)
        y_new += Dn[n - 1] * np.cos(2 * np.pi * n * t_new / T - phi[n - 1])

    plt.figure(figsize=(10, 6))
    plt.plot(t_coords, y_coords, label='Original Data', marker='o')
    plt.plot(t_new, y_fourier, label='Fourier Series Approximation', linestyle='--')
    plt.plot(t_new, y_new, label='Cosine Fourier Series Approximation', linestyle=':')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title('Comparison of Original Data and Fourier Series Approximation')
    plt.legend()
    plt.show()

    return Dn, phi


def write_coefficients_to_excel(file_path, sheet_name, Dn, phi):
    # Load the existing excel file
    book = load_workbook(file_path)

    # Load the sheet into a DataFrame
    df = pd.read_excel(file_path, sheet_name=sheet_name, engine='openpyxl')

    # Find the column index to start writing Dn and phi (one column space from existing content)
    start_col = df.shape[1] + 1

    # Create a new DataFrame to hold the coefficients
    df_coefficients = pd.DataFrame({
        'Dn': Dn,
        'phi': phi
    })

    # Write the coefficients to the Excel sheet
    with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='overlay') as writer:
        df_coefficients.to_excel(writer, sheet_name=sheet_name, index=False, startcol=start_col)

file_path = '/hpc/bsha219/lung/Data/CTEPH/Alfred_Wave_data/digitised_waveforms.xlsx'  # Replace with your file path
sheet_name = 'Alfred9_Post'  # Replace with your sheet name

t_coords, y_coords = read_points_from_excel(file_path, sheet_name)

Y_interpolated = np.interp(np.linspace(min(t_coords), max(t_coords), 100), t_coords, y_coords)

num_coefficients = 20  # Number of terms
A0, An, Bn = compute_fft_coefficients(Y_interpolated, num_coefficients)

Dn, phi = plot_fourier_series(t_coords, y_coords, A0, An, Bn, num_coefficients)
print("A0:", A0)
print("An:", An)
print("Bn:", Bn)
print("Dn:", Dn)
print("phi:", phi)

# Write Dn and phi coefficients to the Excel file
try:
    write_coefficients_to_excel(file_path, sheet_name, Dn, phi)
except Exception as e:
    print(f"Error writing to the Excel file: {e}")
    exit(1)
