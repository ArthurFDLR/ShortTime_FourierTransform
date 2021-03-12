import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

spectrogram_path = './spectrograms/spectrogram.csv'
config_df = pd.read_csv(spectrogram_path, header=None, nrows=1)
config = {}
for entry in config_df.to_numpy()[0]:
    key, value = entry.split(':')
    config[key] = value
print(config)

freq_lim = 3000. #Hz
df = pd.read_csv(spectrogram_path, sep=',',header=None, skiprows=1)
spectrogram = df.to_numpy()[:,0:int(freq_lim/float(config['delta_f']))].transpose()

rcstyle = {'axes.linewidth': 1.0, 'axes.edgecolor': 'black','ytick.minor.size': 5.0}
sns.set(font_scale=2.0)
sns.set_style('ticks', rcstyle)

plt.figure(figsize=(18, 10))
colors = sns.color_palette("Blues", as_cmap=True)
ax = sns.heatmap(spectrogram, cmap=colors, cbar = False)

ax.set_xlabel('Time [s]')
ax.set_xticks(np.linspace(0, spectrogram.shape[1], 10))
ax.set_xticklabels(["{0:.1f}".format(x) for x in np.linspace(0.0, float(config['delta_t']) * spectrogram.shape[1], 10)])

ax.set_ylabel('Frequency [Hz]')
ax.set_yticks(np.linspace(0, 1000, 10))
ax.set_yticklabels([str(int(x)) for x in np.linspace(0.0, float(config['delta_f'])  * freq_lim, 10)])
ax.invert_yaxis()

for side in ['left', 'right', 'bottom', 'top']:
    ax.spines[side].set_visible(True)

plt.tight_layout()
file_name = config['path'].split('/')[-1].split('.')[0]
file_path = './spectrograms/{}_spectrogram.png'.format(file_name)
plt.savefig(file_path)

print("\t- Spectrogram successfully exported: " + file_path)