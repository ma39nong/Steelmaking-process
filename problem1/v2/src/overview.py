import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import *
from interpolation import load_surface
fig=plt.figure(figsize=(17,19))
for k,h in enumerate(HEATS):
 m=load_surface(P/f'data/{h}_model.npz');t=np.linspace(m.time[0],m.time[-1],70);w=np.linspace(m.wave[0],m.wave[-1],120);T,W=np.meshgrid(t,w,indexing='ij');ax=fig.add_subplot(5,2,k+1,projection='3d');ax.plot_surface(T,W,m.evaluate(t,w),cmap='viridis',linewidth=0,rasterized=True);ax.set(title=h,xlabel='Time (s)',ylabel='Wavelength (nm)',zlabel='Intensity');ax.view_init(25,-130)
fig.suptitle('V2 locked B2 + linear interpolation; scientific acceptance incomplete',fontsize=16);fig.tight_layout(rect=(0,0,1,.98));fig.savefig(P/'figures/all_heats_surfaces.png',dpi=150);plt.close(fig)
