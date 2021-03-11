import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

delta_t = 0.0371429
delta_f = 2.69165
freq_lim = int(3000 / 2.69165)

df = pd.read_csv('./spectrograms/chirp.csv', sep=',',header=None)
spectrogram = df.to_numpy()[:,0:freq_lim].transpose()

rcstyle = {'axes.linewidth': 1.0, 'axes.edgecolor': 'black','ytick.minor.size': 5.0}
sns.set(font_scale=2.0)
sns.set_style('ticks', rcstyle)

plt.figure(figsize=(18, 10))
colors = sns.color_palette("Blues", as_cmap=True)
ax = sns.heatmap(spectrogram, cmap=colors, cbar = False)

ax.set_xlabel('Time [s]')
ax.set_xticks(np.linspace(0, spectrogram.shape[1], 10))
ax.set_xticklabels(["{0:.1f}".format(x) for x in np.linspace(0.0, delta_t * spectrogram.shape[1], 10)])

ax.set_ylabel('Frequency [Hz]')
ax.set_yticks(np.linspace(0, 1000, 10))
ax.set_yticklabels([str(int(x)) for x in np.linspace(0.0, delta_f * freq_lim, 10)])
ax.invert_yaxis()

for side in ['left', 'right', 'bottom', 'top']:
    ax.spines[side].set_visible(True)

plt.tight_layout()
plt.savefig('./spectrograms/chirp.png')

print("\t- Spectrogram successfully drawn.")