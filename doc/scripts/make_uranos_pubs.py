# Update URANOS publication list:
# 1. Add new DOIs to publications-uranos.csv 
# 2. Run this script 

# %%
from publib import PubList

# %%
P = PubList()
P.import_csv('publications-uranos.csv').update_citations().sort('date')

# %%
content = """
# Publications using URANOS ({num:.0f})

Citations: **{cites:.0f}** (based on [CrossRef.org](https://www.crossref.org/))

{publist}
""".format(
    num     = P.total_pubs,
    cites   = P.total_cites,
    publist = P.make_list(
        format_str = '- ({year}) {author}  \n**"{title}"**  \n— *{journal}*, doi:[{doi}]({url}), Citations: **{cited}**  \n')
)[1:]

# %%
with open('../PUBLICATIONS.md', 'w', encoding="utf-8") as fh:
    fh.write(content)
    
# %%
from pandas import DataFrame
history = DataFrame([dict(
    year = today.year,
    pubs = P.total_pubs)])

# %%
plotdata = P.data[['year','authors','cites']]
plotdata['firstauthor'] = [x[0] for x in plotdata.authors]
plotdata['label'] = ['{1}\n{0} et al.'.format(x[0], y) for x,y in zip(plotdata.authors, plotdata.year)]
plotdata = plotdata[::-1]
plotdata

# %%
from spiral_strip_library_v04 import *
duration = [ 7,  5,  9, 23,  6,  6,  4]
deaths = [ 85., 28., 22.5, 9.846, 7., 3.95, 3. ]

color1='#e82b33'
color2='#31aae0'
segment_color=[]
for il in range(len(duration)):
	segment_color.append(color1)
	segment_color.append(color2)

segment_color='C0'

width = []
segment_length = 150

draw_spiral_strip(image_filename="example06.pdf", r0=0, space=0, fi0_deg=0,
    number_of_segments = len(plotdata), width=plotdata.cites.values+1, segment_length=segment_length,
    maximum_length_of_subsegment=10, segment_color = segment_color,
    values=None, colormap_name=None,
    labels1_text=plotdata.label.values, labels1_fontsize=5, labels1_color="dimgray", labels1_pad=35,
    antialiased = True, gap_between_consecutive_segments = 15,
    dpi=300)

# %%
from corny.figures import Figure
c_pubs = 'C0'
c_cite = '#000000'
with Figure(size=(8,3), save='out/pubs_per_year.png') as ax:
    ax.bar(history[:-1].year, history[:-1].pubs, color=c_pubs, alpha=0.5)
    ax.bar(history[-1:].year, history[-1:].pubs, color=c_pubs, alpha=0.5, hatch='//', edgecolor='white')
    #ax.scatter(today.year, spy_now.pubs.values[-1] + PP.total_pubs, color='#BBBBBB')
    ax.set_ylabel('Publications')
    ax.spines['left'].set_color(c_pubs)
    ax.spines['top'].set_color('w')
    ax.yaxis.label.set_color(c_pubs)
    ax.tick_params(axis='y', colors=c_pubs)
