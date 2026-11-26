"""Download publicly offered author/publisher texts; verify PDF signatures. No paywall bypass."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests,json,hashlib
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parent
(ROOT/'fulltext').mkdir(exist_ok=True)
(ROOT/'metadata').mkdir(exist_ok=True)
sources={
'R03':['https://arxiv.org/pdf/2201.02812'],
'R09':['https://www.oceanoptics.com/wp-content/uploads/2024/12/USB4000-Operating-Instructions.pdf'],
'R01_cran':['https://cran.r-project.org/web/packages/JOPS/JOPS.pdf'],
'S06':['https://andrewd.ces.clemson.edu/courses/cpsc881/papers/Sch11_whatIsSG.pdf'],
'R01':['https://www.ma.hw.ac.uk/~iain/research/GLAM/CSDA_2006/CSDA.pdf'],
'R04':['https://www.nature.com/articles/s41598-023-41751-9.pdf','https://pmc.ncbi.nlm.nih.gov/articles/PMC10475112/'],
'R05':['https://opg.optica.org/josaa/fulltext.cfm?uri=josaa-31-7-1445&id=288946','https://opg.optica.org/josaa/viewmedia.cfm?uri=josaa-31-7-1445&seq=0'],
'R06':['https://pubs.acs.org/doi/pdf/10.1021/ac60214a047'],
'R07':['https://journals.sagepub.com/doi/full/10.1177/03019233241273473'],
'R08':['https://onlinelibrary.wiley.com/doi/pdf/10.1155/2018/8682725','https://downloads.hindawi.com/journals/complexity/2018/8682725.pdf'],
'R10':['https://eeweb.engineering.nyu.edu/iselesni/pubs/BEADS_2014.pdf','https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/49974/versions/3/previews/BEADS_2014.pdf'],
'R11':['https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2Faoms%2F1177703732','https://projecteuclid.org/euclid.aoms/1177703732'],
'R12':['https://journals.sagepub.com/doi/full/10.1177/03019233251382772']}
def get(rid,urls):
 logs=[]
 for k,u in enumerate(urls):
  try:
   try:r=requests.get(u,timeout=35)
   except requests.exceptions.SSLError:
    r=requests.get(u,timeout=35,verify=False);logs.append({'url':u,'note':'server certificate chain unverified; content hash recorded'})
   r.raise_for_status();pdf=r.content.startswith(b'%PDF');ext='pdf' if pdf else 'html'
   name=rid if rid in ['R03','R09','R01_cran','S06'] and k==0 else f'{rid}_{k}'
   path=ROOT/'fulltext'/f'{name}.{ext}';path.write_bytes(r.content)
   if pdf:
    from pypdf import PdfReader
    path.with_suffix('.txt').write_text('\n'.join(f'=== PDF PAGE {i+1} ===\n'+(page.extract_text() or '') for i,page in enumerate(PdfReader(path).pages)))
   if not pdf:(ROOT/'fulltext'/f'{rid}_{k}.txt').write_text(BeautifulSoup(r.text,'html.parser').get_text('\n',strip=True))
   logs.append({'url':u,'status':r.status_code,'file':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(r.content).hexdigest(),'bytes':len(r.content),'pdf':pdf})
   print(rid,ext,len(r.content),flush=True)
  except Exception as e:logs.append({'url':u,'error':str(e)});print(rid,str(e)[:180],flush=True)
 (ROOT/'metadata'/f'{rid}_acquisition.json').write_text(json.dumps(logs,ensure_ascii=False,indent=2))
with ThreadPoolExecutor(6) as ex:list(ex.map(lambda x:get(*x),sources.items()))
# Author publication list: derive exact R02 PDF link rather than inventing its filename.
u='https://numana.uni-koeln.de/arbeitsgruppe/lehrstuhlinhaberin/prof-dr-angela-kunoth/publications-and-reports'
try:
 r=requests.get(u,timeout=45);s=BeautifulSoup(r.text,'html.parser')
 for a in s.find_all('a',href=True):
  if 'Robust Regression' in a.get_text() or 'robust' in a['href'].lower():
   from urllib.parse import urljoin
   print('R02 found',a.get_text(),a['href']);get('R02',[urljoin(u,a['href'])])
except Exception as e:print(e)
