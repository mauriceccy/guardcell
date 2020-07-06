import scobra
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#     Loading CSV
scan = pd.read_csv("GC_sucrose_scan_standard.csv",header=None,index_col=0)

#     Extracting Sucrose_tx_Day data
sucrose_tx_flux = scan.loc['Sucrose_tx_Day']/11.0*1000

#     1. Photon
photon_tx_flux = scan.loc['Photon_tx_Open'] + scan.loc['Photon_tx_Day']

#     2. Rubisco carboxylase
rubisco_day_flux = scan.loc['RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Day']/11.0

#     3. OPPP
g6pdh_day_flux = (scan.loc['GLU6PDEHYDROG_RXN_c_Day'] + scan.loc['GLU6PDEHYDROG_RXN_p_Day']) / 11.0

#     4. Lower glycolysis
pk_day_flux = (scan.loc['PEPDEPHOS_RXN_c_Day'] + scan.loc['PEPDEPHOS_RXN_p_Day']) / 11.0

#     5. TCA
complex2_day_flux = scan.loc['SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN_mi_Day'] / 11.0

#     6. Starch accumulation
starch_dayclose_flux = scan.loc['Starch_DayClose_storage']
starch_closenight_flux = scan.loc['Starch_CloseNight_storage']
starch_nightopen_flux = scan.loc['Starch_NightOpen_storage']
starch_openday_flux = scan.loc['Starch_OpenDay_storage']


#     Creating the figure
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(3, 2, figsize=(8,8))
fig.subplots_adjust(top=2.5, right=2, hspace=.35, wspace=.3)

# fonts
axes_fontsize = 20
index_fontsize = 24
xlab_fontsize = 18
ylab_fontsize = 18
legend_fontsize = 18

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
sucroseunit = r'$(\mu mol\; g^{-1} DW\;h^{-1})$'
dayunit = r'$(mmol\; g^{-1} DW\;day^{-1})$'
starchunit = r'$(mmol\; g^{-1} DW)$'

# 1. Photon
ax[0,0].scatter(sucrose_tx_flux, photon_tx_flux,s=10,c=colour1, marker=shape1)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='Total photon influx')]
ax[0,0].legend(handles=legend_elements, prop={'size': legend_fontsize})
ax[0,0].set_xlabel('Sucrose uptake in Day phase ' + sucroseunit, fontsize=axes_fontsize)
ax[0,0].set_ylabel('Total photon influx \n' + dayunit, fontsize=axes_fontsize)
ax[0,0].tick_params(axis="x", labelsize=xlab_fontsize)
ax[0,0].tick_params(axis="y", labelsize=ylab_fontsize)
ax[0,0].text(-0.2, 1.1, '(a)', transform=ax[0,0].transAxes, size=index_fontsize, weight='bold')

# 2. Rubisco carboxylase
ax[0,1].scatter(sucrose_tx_flux, rubisco_day_flux,s=10,c=colour1, marker=shape1)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='Rubisco carboxylase flux')]
ax[0,1].legend(handles=legend_elements, prop={'size': legend_fontsize})
ax[0,1].set_xlabel('Sucrose uptake in Day phase ' + sucroseunit, fontsize=axes_fontsize)
ax[0,1].set_ylabel('Rubisco carboxylase flux in Day phase \n' + unit, fontsize=axes_fontsize)
ax[0,1].tick_params(axis="x", labelsize=xlab_fontsize)
ax[0,1].tick_params(axis="y", labelsize=ylab_fontsize)
ax[0,1].text(-0.2, 1.1, '(b)', transform=ax[0,1].transAxes, size=index_fontsize, weight='bold')

# 3. OPPP
ax[1,0].scatter(sucrose_tx_flux, g6pdh_day_flux,s=10,c=colour1, marker=shape1)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='OPPP flux')]
ax[1,0].legend(handles=legend_elements, prop={'size': legend_fontsize})
ax[1,0].set_xlabel('Sucrose uptake in Day phase ' + sucroseunit, fontsize=axes_fontsize)
ax[1,0].set_ylabel('OPPP flux in Day phase \n' + unit, fontsize=axes_fontsize)
ax[1,0].tick_params(axis="x", labelsize=xlab_fontsize)
ax[1,0].tick_params(axis="y", labelsize=ylab_fontsize)
ax[1,0].text(-0.2, 1.1, '(c)', transform=ax[1,0].transAxes, size=index_fontsize, weight='bold')


# 4. lower glycolysis
ax[1,1].scatter(sucrose_tx_flux, pk_day_flux,s=10,c=colour1, marker=shape1)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='Lower glycolysis flux')]
ax[1,1].legend(handles=legend_elements, prop={'size': legend_fontsize})
ax[1,1].set_xlabel('Sucrose uptake in Day phase ' + sucroseunit, fontsize=axes_fontsize)
ax[1,1].set_ylabel('Lower glycolysis flux in Day phase\n' + unit, fontsize=axes_fontsize)
ax[1,1].tick_params(axis="x", labelsize=xlab_fontsize)
ax[1,1].tick_params(axis="y", labelsize=ylab_fontsize)
ax[1,1].text(-0.2, 1.1, '(d)', transform=ax[1,1].transAxes, size=index_fontsize, weight='bold')

# 5. TCA
ax[2,0].scatter(sucrose_tx_flux, complex2_day_flux,s=10,c=colour1, marker=shape1)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='TCA cycle flux')]
ax[2,0].legend(handles=legend_elements, prop={'size': legend_fontsize})
ax[2,0].set_xlabel('Sucrose uptake in Day phase ' + sucroseunit, fontsize=axes_fontsize)
ax[2,0].set_ylabel('TCA cycle flux in Day phase\n' + unit, fontsize=axes_fontsize)
ax[2,0].tick_params(axis="x", labelsize=xlab_fontsize)
ax[2,0].tick_params(axis="y", labelsize=ylab_fontsize)
ax[2,0].text(-0.2, 1.1, '(e)', transform=ax[2,0].transAxes, size=index_fontsize, weight='bold')

# 6. Starch accumulation
ax[2,1].scatter(sucrose_tx_flux,starch_dayclose_flux,s=10,c=colour1, marker=shape1)
ax[2,1].scatter(sucrose_tx_flux,starch_closenight_flux,s=10,c=colour2, marker=shape2)
ax[2,1].scatter(sucrose_tx_flux,starch_nightopen_flux,s=10,c=colour3, marker=shape3)
ax[2,1].scatter(sucrose_tx_flux,starch_openday_flux,s=10,c=colour4, marker=shape4)
legend_elements = [Line2D([0], [0], color=colour1, marker=shape1, label='Day-Close'),
                   Line2D([0], [0], color=colour2, marker=shape2, label='Close-Night'),
                   Line2D([0], [0], color=colour3, marker=shape3, label='Night-Open'),
                   Line2D([0], [0], color=colour4, marker=shape4, label='Open-Day')]
ax[2,1].legend(handles=legend_elements, prop={'size': legend_fontsize})
ax[2,1].set_xlabel('Sucrose uptake in Day phase ' + sucroseunit, fontsize=axes_fontsize)
ax[2,1].set_ylabel('Starch accumulation \n ' + starchunit, fontsize=axes_fontsize)
ax[2,1].tick_params(axis="x", labelsize=xlab_fontsize)
ax[2,1].tick_params(axis="y", labelsize=ylab_fontsize)
ax[2,1].text(-0.2, 1.1, '(f)', transform=ax[2,1].transAxes, size=index_fontsize, weight='bold')

fig.savefig('Figure_3.png',bbox_inches = 'tight' ,dpi=300)
#plt.show()

