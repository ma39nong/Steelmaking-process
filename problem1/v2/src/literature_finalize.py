from common import *
import pandas as pd
R=P/'references';d=pd.read_csv(R/'literature_index_v2.csv',keep_default_na=False);bib=[]
for i,row in d.iterrows():
 id=row.id;mp=R/f'{id}_metadata.json';m=json.loads(mp.read_text()) if mp.exists() else {}
 journal=(m.get('container-title') or ['Journal of Machine Learning Research'])[0];volume=m.get('volume','11');issue=m.get('issue','80' if id=='V14' else '');pages=m.get('page','2287–2322' if id=='V14' else '')
 d.loc[i,'journal']=journal;d.loc[i,'volume']=volume;d.loc[i,'issue']=issue;d.loc[i,'pages_published']=pages
 note=R/'notes'/f'{id}.md';s=note.read_text();line=f'书目信息：{journal}, {volume}({issue}): {pages}。'
 if line not in s:s=s.replace('\n\n作者：','\n\n'+line+'\n\n作者：');note.write_text(s)
 fields=dict(title=row.title,author=row.authors,year=int(row.year),journal=journal,volume=volume,number=issue,pages=pages,url=row.url)
 if row.doi:fields['doi']=row.doi
 bib.append('@article{'+id+',\n'+',\n'.join(' '+k+'={'+str(v)+'}' for k,v in fields.items() if v!='')+'\n}')
d.to_csv(R/'literature_index_v2.csv',index=False);(R/'references_v2.bib').write_text('\n\n'.join(bib)+'\n')
