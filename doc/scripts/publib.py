# Collect publication details and citations based on a list of DOIs and CrossRef.org
# v1.0 by M. Schrön (2022)

import pandas
from habanero import Crossref, counts
from datetime import date

"""
CrossRef columns
'indexed', 'reference-count', 'publisher', 'issue',
'license', 'content-domain', 'short-container-title',
'DOI', 'type', 'created', 'page', 'update-policy',
'source', 'is-referenced-by-count', 'title', 'prefix',
'volume', 'author', 'member', 'published-online',
'reference', 'container-title', 'original-title',
'language', 'link', 'deposited', 'score', 'resource',
'subtitle', 'editor', 'short-title', 'issued',
'journal-issue', 'URL', 'relation', 'ISSN',
'issn-type', 'subject', 'published'    
"""


class CRPub:
    def __init__(self, CRitem, me='Einstein'):
        self.data = CRitem['message']
        
    def get_date(self):
        date_items = self.data['published']['date-parts'][0]
        if len(date_items) < 3 and 'published-online' in self.data:
            date_items = self.data['published-online']['date-parts'][0]
        if len(date_items) == 2:
            date_items = date_items[0], date_items[1], 1
        return(date(*date_items))

    def get_year(self):
        date_items = self.data['published']['date-parts'][0]
        year = int(date_items[0])
        return(year)

    def get_title(self):
        title = self.data['title']
        if isinstance(title, list):
            title = title[0]
        return(title)
        
    def get_authors(self):
        author_dict_list = self.data['author']
        authors = ['%s' % author['family'] for author in author_dict_list]
        return(authors)

    def get_cites(self):
        cites = counts.citation_count(doi=self.get_doi())
        return(cites)

    def get_doi(self):
        doi = self.data['DOI']
        return(doi)

    def get_url(self):
        url = self.data['URL']
        return(url)
        
    def get_journal(self):
        if 'subtype' in self.data and self.data['subtype']=='preprint':
            journal='[preprint]'
        else:
            if 'short-container-title' in self.data:
                journal = self.data['short-container-title']
            if not journal:
                journal = self.data['container-title']
        if isinstance(journal, list) and len(journal)>0:
            journal = journal[0]
        return(journal)
        
    
class PubList:
    
    def __init__(self, me='Einstein'):
        self.me = me
        self.CR = Crossref()
        
    def import_csv(self, file):
        self.data = pandas.read_csv(file, skipinitialspace=True)
        self._CRdata = self.CR.works(ids=self.data.doi.values)
        self.Pubs = [CRPub(CRitem, me=self.me) for CRitem in self._CRdata]
        self.data['date']    = [pub.get_date()    for pub in self.Pubs]
        self.data['title']   = [pub.get_title()   for pub in self.Pubs]
        self.data['authors'] = [pub.get_authors() for pub in self.Pubs]
        self.data['year']    = [pub.get_year()    for pub in self.Pubs]
        self.data['doi']     = [pub.get_doi()     for pub in self.Pubs]
        self.data['url']     = [pub.get_url()     for pub in self.Pubs]
        self.data['journal'] = [pub.get_journal()   for pub in self.Pubs]
        if 'team_num' in self.data:
            self.data['team_num']      = self.data['team_num'].fillna(0)
        if 'team_notes' in self.data:
            self.data['team_notes']    = self.data['team_notes'].fillna('')
        if 'corresponding' in self.data:
            self.data['corresponding'] = self.data['corresponding'].fillna('')
        
        self.total_pubs    = len(self.data)
        self.date_earliest = self.data.date.min() #.strftime('%Y %b')
        self.date_latest   = self.data.date.max() #.strftime('%Y %b')
        return(self)
        
    def update_citations(self):
        self.data['cites'] = [pub.get_cites() for pub in self.Pubs]
        self.total_cites = self.data.cites.sum()
        return(self)

    def sort(self, column='date', ascending=False):
        self.data = self.data.sort_values(column, ascending=ascending)
        return(self)
        
    def make_list(self,
        format_str = '- ({year}) {author}  \n"{title}"  \n— *{journal}*, [doi:{doi}]({url}), Citations: **{cited}**  \n',
        authors_kw = dict(max_authors=10, highlight_me='bold')):
        
        r = ''
        for i, row in self.data.iterrows():
            r += format_str.format(
               title   = row['title'],
               year    = row['year'],
               author  = self.format_authors(row['authors'],
                    corresponding=row['corresponding'],
                    team_num=row['team_num'],
                    team_notes=row['team_notes'],
                    **authors_kw),
               journal = row['journal'],
               cited   = row['cites'],
               doi     = row['doi'],
               url     = row['url'] )
        return(r)
    
    def format_authors(self, authors, join=', ', corresponding=False, 
        max_authors=None, highlight_me='italic', export='md',
        team_num=None, team_notes=None):

        if not team_num:
            team_num = 0
        r = authors
        if not max_authors is None:
            i_me = authors.index(self.me) if self.me in authors else -1
            i_last = len(authors)-1
            is_last = True if i_me == i_last else False
            xauthors = []
            i = 0
            includes_me = False
            
            for author in authors:
                if author == self.me:
                    includes_me = True
                if i < team_num:
                    author += '\*'
            
                if i==0 or i==i_me or i==i_last:
                    xauthors.append(author)
                elif i < i_last and includes_me and i < max_authors-1:
                    xauthors.append(author)
                elif i < i_last and not includes_me and i < max_authors-1 and is_last:
                    xauthors.append(author)
                elif i < i_last and not includes_me and i < max_authors-2 and not is_last:
                    xauthors.append(author)
                else:
                    if xauthors[-1][0] == '.':
                        xauthors[-1] += '.'
                    else:
                        xauthors.append('.')
                i +=1
            r = xauthors

        MEf = self.me
        if corresponding:
            r = list(map(lambda x: x.replace(MEf, MEf+'\*\*'), r))
                
        if highlight_me:
            if highlight_me == 'italic':
                if export=='md':   MEf = '*%s*' % MEf
                if export=='html': MEf = '<i>%s</i>' % MEf
                if export=='tex':  MEf = r'\emph{%s}' % MEf
            elif highlight_me == 'bold':
                if export=='md':   MEf = '**%s**' % MEf
                if export=='html': MEf = '<b>%s</b>' % MEf
                if export=='tex':  MEf = r'\textbf{%s}' % MEf
            elif highlight_me == 'underline':
                if export=='md':   MEf = '_%s_' % MEf
                if export=='html': MEf = '<u>%s</u>' % MEf
                if export=='tex':  MEf = r'\underline{%s}' % MEf

            r = list(map(lambda x: x.replace(self.me, MEf), r))

        r = join.join(r)
        
        if team_num > 0:
            r += ' *(\*%s)*' % team_notes
        if corresponding:
            r += ' *(\*\*corresponding)*'
        
        return(r)
