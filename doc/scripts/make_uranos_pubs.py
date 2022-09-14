# Update URANOS publication list:
# 1. Add new DOIs to publications-uranos.csv 
# 2. Run this script 

# %%
from publib import PubList
from datetime import date

# %%
P = PubList()
P.import_csv('publications-uranos.csv').update_citations().sort('date')

# %%
import compress_pickle as pickle
with open('pubdata-cache.bz2', 'wb') as file:
    pickle.dump(P, file)
with open('pubdata-cache.bz2', 'rb') as file:
    P = pickle.load(file)

# %%
content = """
# Publications using URANOS ({num:.0f})

Citations: **{cites:.0f}** (based on [CrossRef.org](https://www.crossref.org/))

*Figure. Left: Publications and their total citations (bar width). Right: Cumulative publications over the years. (Spiral package credits: [G. Skok, 2022](https://doi.org/10.3390/app12136609))*
![Publications and citations per year](pubplot.png)

## Details 
{publist}
""".format(
    num     = P.total_pubs,
    cites   = P.total_cites,
    publist = P.make_list(
        format_str = '- `{year}` {author}  \n**"{title}"**  \n— *{journal}*, doi:[{doi}]({url}), Citations: **{cited}**  \n')
)[1:]

# %%
with open('../PUBLICATIONS.md', 'w', encoding="utf-8") as fh:
    fh.write(content)
    
# %%
import numpy as np
from pandas import DataFrame, to_datetime
plotdata = P.data[['year','date','authors','cites']]
plotdata['firstauthor'] = [x[0] for x in plotdata.authors]
plotdata['label'] = ['{1}\n{0} et al.'.format(x[0], y) for x,y in zip(plotdata.authors, plotdata.year)]
plotdata = plotdata[::-1]
plotdata['t'] = to_datetime(plotdata.date)
plotdata['num'] = 1
plotdata['num'] = np.cumsum(plotdata['num'])
plotdata

# %%
import matplotlib.pyplot as plt
# Corny package from: https://git.ufz.de/CRNS/cornish_pasdy
from corny.figures import Figure
# Sprial package from: https://github.com/skokg/Spiral-strip
from spiral_strip_library_v04 import *

radius          = 0
space           = 0
rotation        = 110 # try to keep the latest paper on top
segment_len     = 130
segment_len_max = 10
segment_gap     = 10
label_padding   = 25
label_fontsize  = 5

color_data   = '#1f77b4'
color_labels = '#BBBBBB'

#themes = dict(#
#    dark  = dict(style='default',         nametag='-light', fgcolor='black', datacolor='#1f77b4'),
#    light = dict(style='dark_background', nametag='-dark',  fgcolor='white', datacolor='#1f77b4'))

#58a6ff
#for theme in themes:
#    plt.style.use(themes[theme]['style'])
with Figure(layout=(1,2), size=(8,4),
    save='../pubplot.png') as axes:
    
    ax = axes[0]
    draw_spiral_strip(r0=radius, space=space, fi0_deg=rotation,
        number_of_segments=len(plotdata), width=plotdata.cites.values+1,
        segment_length=segment_len, maximum_length_of_subsegment=segment_len_max,
        segment_color=color_data, values=None, colormap_name=None, labels1_color=color_labels,
        labels1_text=plotdata.label.values, labels1_fontsize=label_fontsize,
        labels1_pad=label_padding, antialiased=True,
        gap_between_consecutive_segments = segment_gap,
        ax=ax)
    ax.autoscale(enable=True, axis='x', tight=True)
    
    ax = axes[1]
    ax.set_title('Publications (total)', fontsize=7, color=color_labels)
    ax.plot(plotdata.t, plotdata.num, color=color_data, drawstyle='steps-post')
    ax.set_xlim(date(2015,1,1), date(2022,12,31))
    ax.set_yticks(np.arange(0,len(plotdata),5))
    ax.xaxis.set_tick_params(labelsize=6)
    ax.yaxis.set_tick_params(labelsize=6)
    ax.grid(color=color_labels, ls=':', axis='y', alpha=1)
    ax.spines.top.set(visible=False)
    ax.spines.left.set(visible=False)
    ax.spines.right.set(visible=False)
    ax.tick_params(colors=color_labels, which='both', axis='y', tick1On=False)
    ax.tick_params(colors=color_labels, which='both', axis='x')
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["bottom"].set_color(color_labels)
    ax.plot(1, 0, ">", color=color_labels, transform=ax.get_yaxis_transform(), markersize=3, clip_on=False)
