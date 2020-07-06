import scobra
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#     Loading CSV
scan = pd.read_csv("GC_ATP_scan_standard.csv",header=None,index_col=0)

#     Extracting ATPase_tx_data
atpase_tx_flux = scan.loc['ATPase_tx_Day']/11.0

#     Starch accumulation
starch_dayclose_flux = scan.loc['Starch_DayClose_storage']
starch_closenight_flux = scan.loc['Starch_CloseNight_storage']
starch_nightopen_flux = scan.loc['Starch_NightOpen_storage']
starch_openday_flux = scan.loc['Starch_OpenDay_storage']


#     Creating the figure
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

#fig, ax = plt.subplots(2, 2, figsize=(8,8))
#fig.subplots_adjust(top=2.5, right=2, hspace=.35, wspace=.3)

# fonts
axes_fontsize = 12
index_fontsize = 12
xlab_fontsize = 12
ylab_fontsize = 12
legend_fontsize = 12

# colour codes can be found at https://matplotlib.org/examples/color/named_colors.html
colour1='blue'
colour2='maroon'
colour3 = 'orangered'
colour4 = 'darkcyan'
# shape codes can be found at https://matplotlib.org/3.1.0/api/markers_api.html#module-matplotlib.markers
shape1='o'
shape2='^'
shape3='x'
shape4='P'
# unit
unit = r'$(mmol\; g^{-1} DW\;h^{-1})$'
atpunit = r'$(mmol\; g^{-1} DW\;h^{-1})$'
dayunit = r'$(mmol\; g^{-1} DW\;day^{-1})$'
starchunit = r'$(mmol\; g^{-1} DW)$'

# Starch accumulation
plt.scatter(atpase_tx_flux,starch_dayclose_flux,s=10,c=colour1, marker=shape1)
plt.scatter(atpase_tx_flux,starch_closenight_flux,s=10,c=colour2, marker=shape2)
plt.scatter(atpase_tx_flux,starch_nightopen_flux,s=10,c=colour3, marker=shape3)
plt.scatter(atpase_tx_flux,starch_openday_flux,s=10,c=colour4, marker=shape4)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='Day-Close'),
                   Line2D([0], [0], color=colour2, marker=shape2, label='Close-Night'),
                   Line2D([0], [0], color=colour3, marker=shape3, label='Night-Open'),
                   Line2D([0], [0], color=colour4, marker=shape4, label='Open-Day')]
plt.legend(handles=legend_elements, prop={'size': legend_fontsize})
plt.xlabel('ATP maintenance ' + atpunit, fontsize=axes_fontsize)
plt.ylabel('Starch accumulation \n ' + starchunit, fontsize=axes_fontsize)
plt.tick_params(axis="x", labelsize=xlab_fontsize)
plt.tick_params(axis="y", labelsize=ylab_fontsize)
#ax[0,0].text(-0.2, 1.1, '(f)', transform=ax[2,1].transAxes, size=index_fontsize, weight='bold')

plt.savefig('Figure_S1.png',bbox_inches = 'tight' ,dpi=300)
#plt.show()
