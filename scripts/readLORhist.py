import numpy as np
import matplotlib.pyplot as plt

file_path = 'LORhist.dat'

data_mm = np.memmap(file_path, dtype='float32', mode='r')

step = 10
sampled_data = data_mm[::step]

plt.figure(figsize=(12, 6))
plt.plot(sampled_data)
plt.title(f'Histograma LOR (Muestreado 1:{step})')
plt.xlabel(f'Bin Index (x{step})')
plt.ylabel('Cuentas')
plt.grid(True)
plt.show()