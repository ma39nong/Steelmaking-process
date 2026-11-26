from pathlib import Path
import requests,json,hashlib
from pypdf import PdfReader
P=Path(__file__).resolve().parents[1]/'references'
items=[('V13','10.1145/1970392.1970395','https://www.columbia.edu/~jw2966/papers/CLMW11-JACM.pdf'),('V14','10.5555/1756006.1859931','https://jmlr.org/papers/volume11/mazumder10a/mazumder10a.pdf'),('V15','10.1137/0717021','https://epubs.siam.org/doi/pdf/10.1137/0717021'),('V16','10.1137/0905021','https://epubs.siam.org/doi/pdf/10.1137/0905021'),('V17','10.1093/imaiai/iau006','https://arxiv.org/pdf/1209.3672')]
rows=[]
for id,doi,url in items:
 row=dict(id=id,doi=doi,url=url)
 try:
  r=requests.get('https://api.crossref.org/works/'+doi,timeout=40);row['crossref_status']=r.status_code
  if r.ok:(P/f'{id}_metadata.json').write_text(json.dumps(r.json()['message'],ensure_ascii=False,indent=2))
 except Exception as e:row['metadata_error']=str(e)
 try:
  r=requests.get(url,timeout=50);row['http_status']=r.status_code
  if r.ok and r.content.startswith(b'%PDF'):
   f=P/'fulltext'/f'{id}.pdf';f.write_bytes(r.content);doc=PdfReader(f);row.update(status='公开全文已获取',pages=len(doc.pages),sha256=hashlib.sha256(r.content).hexdigest());(P/'fulltext'/f'{id}.txt').write_text('\n'.join(f'\n=== PDF PAGE {i+1} ===\n'+(p.extract_text() or '') for i,p in enumerate(doc.pages)))
  else:row['status']='未获取合法全文；不声称阅读全文'
 except Exception as e:row.update(status='获取失败',error=str(e))
 rows.append(row);print(row,flush=True)
(P/'retrieval.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
