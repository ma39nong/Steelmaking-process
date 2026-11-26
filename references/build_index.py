from pathlib import Path
import json,csv,hashlib,html
R=Path(__file__).resolve().parent
status={1:('未取得全文；出版社预览已读',''),2:('未取得全文；PubMed摘要及作者目录已核实',''),3:('已下载并阅读作者预印本；与最终版作者列表不同','fulltext/R03.pdf'),4:('已下载并阅读出版社OA全文','fulltext/R04_0.pdf'),5:('未取得全文；出版社摘要与作者主页已核实',''),6:('原文付费未取得；补充阅读S06全文',''),7:('未取得全文；仅书目与摘要核实',''),8:('已阅读出版社HTML正文；PDF下载被拒','fulltext/R08_web_extracted.txt'),9:('已下载并阅读厂家官方手册','fulltext/R09.pdf'),10:('已下载并阅读作者公开预印本','fulltext/R10_0.pdf'),11:('公开入口被访问验证拦截；原文未取得',''),12:('未取得全文；仅出版社书目与摘要核实','')}
roles={1:'张量积加权P样条的主模型来源候选；全文缺口由JOPS作者文档与独立推导补充',2:'局部支撑与异常影响分析；不冒充原算法复现',3:'低秩稀疏/TV对照的思想来源；二维改造',4:'炼钢背景/数据条件/BEADS独立实验',5:'整条异常光谱分组处理；PCA后MCD是本题改造',6:'传统局部多项式平滑基线',7:'后续烟气预测储备，本问仅背景',8:'非平稳与多尺度背景',9:'量程、校准、仪器噪声与饱和假设',10:'BEADS目标与求解的直接算法依据',11:'Huber理论来源候选；实现公式由官方SciPy文档核实',12:'后续问题储备；不复现网络'}
rows=[];bib=[]
for k in range(1,13):
 rid=f'R{k:02d}';p=R/f'metadata/{rid}_crossref.json'
 if p.exists():
  d=json.loads(p.read_text());authors=' and '.join(x.get('given','')+' '+x['family'] for x in d.get('author',[]));title=d['title'][0];year=d['published']['date-parts'][0][0];doi=d['DOI'];journal=html.unescape(d.get('container-title',[''])[0]);vol=d.get('volume','');issue=d.get('issue','');pages=d.get('page',d.get('article-number',''))
  if k==5:pages='1445--1452'
  if k==4:pages='14409'
  if k==8:pages='8682725'
 else:
  authors='Ocean Optics, Inc.';title='USB4000 Fiber Optic Spectrometer Installation and Operation Manual';year=2016;doi='';journal='Document 211-00000-000-02-201604';vol=issue=pages=''
 if k==7:year=2026
 url='https://doi.org/'+doi if doi else 'https://www.oceanoptics.com/wp-content/uploads/2024/12/USB4000-Operating-Instructions.pdf'
 state,file=status[k]
 if file and not (R/file).exists():state+='；当前本机缓存不存在，需重新获取'
 sha=hashlib.sha256((R/file).read_bytes()).hexdigest() if file and (R/file).exists() else ''
 row=dict(id=rid,title=title,authors=authors,year=year,journal=journal,volume=vol,issue=issue,pages=pages,doi=doi,url=url,access_status=state,local_file=file,sha256=sha,role=roles[k],access_date='2026-11-26');rows.append(row)
 entry='@manual' if k==9 else '@article'
 fields={'title':title,'author':authors,'year':str(year),'journal':journal,'volume':vol,'number':issue,'pages':str(pages).replace('--','-').replace('-','--'),'doi':doi,'url':url,'note':state+('; first published online 2024-08-22, print 2026-03' if k==7 else '')+'; accessed 2026-11-26'}
 bib.append(entry+'{'+rid+',\n'+',\n'.join('  '+f+' = {'+v+'}' for f,v in fields.items() if v)+'\n}\n')
with (R/'literature_index.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(R/'references.bib').write_text('\n'.join(bib)+(R/'supplementary.bib').read_text())
(R/'metadata/index.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
