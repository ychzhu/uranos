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
    
    
