# Update URANOS publication list:
# 1. Add new DOIs to publications-uranos.csv 
# 2. Run this script 

# %%
from publib import PubList

# %%
P = PubList()
P.import_csv('publications-uranos.csv').update_citations().sort('date')

# %%
with open('../PUBLICATIONS.md', 'w', encoding="utf-8") as fh:
    fh.write('# Publications using URANOS ({num:.0f})\n\n'.format(
        num   = P.total_pubs ))
    fh.write(P.make_list())
    
