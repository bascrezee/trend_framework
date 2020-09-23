import numpy as np
import iris
import iris.plot
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import os

def trend_mapplot(indata, args):
    plt.clf()
    fig = plt.figure(dpi=200)
    fig.add_subplot(projection=ccrs.PlateCarree())

    plottitle =  args.datasetname

    #plotkwargs = _parse_cmap(plotkwargs) Check
    extend = 'both'

    # Prepare cmap
    nbins = 11 # has to be uneven
    cmapname = 'BrBG'
    vmin = -0.001
    vmax = 0.001
    cmap = matplotlib.cm.get_cmap(cmapname)
    color_list = cmap(np.linspace(0, 1, nbins))
    cmap_name = cmap.name + str(nbins)
    cmap = cmap.from_list(cmap_name, color_list, nbins)

    levels = MaxNLocator(nbins=nbins, symmetric=True).tick_values(
        vmin, vmax)
    # Remove zero from levels
    levels = np.delete(levels, len(levels) / 2)

    color_list = list(color_list)

    cmap, norm = from_levels_and_colors(levels, color_list,
                                        extend=extend)
    cmap.set_bad("grey", 0.1)

    pcols = indata.plot.pcolormesh(transform=ccrs.PlateCarree() , cmap=cmap, norm=norm, rasterized=True, add_colorbar=False)
    # Take out small grid lines like this
    pcols.set_edgecolor('face')
    plt.gca().coastlines()

    # Colorbar
    colorbar = plt.colorbar(pcols, orientation='horizontal', extend='both')
    colorbar.set_label(indata.attrs['units'])
    colorbar.ax.tick_params(labelsize=8)

    # Add timeperiod to plot title
    # Get timeperiod from args
    a = args.startdate.split('-')[0]
    b = args.enddate.split('-')[0]
    timeperiod = f"{a}-{b}"
    plt.title(f"{plottitle} {timeperiod}")
    savename = os.path.join(args.outputdir, f"theilsenmk_{plottitle}_{a}_{b}.png")
    print(f"Saving {savename}")
    fig.savefig(savename, bbox_inches='tight')
    plt.close(fig)
