"""Prepare figure-linked manuscript and verify targeted visual redesign."""
from pathlib import Path
import re,json,hashlib,subprocess
import pandas as pd,numpy as np
from PIL import Image,ImageDraw
O=Path(__file__).resolve().parents[1];F=O.parents[1];R=F/'reports';checks={}
def verify(name,yes):
 assert yes,name
 checks[name]=True
caps={
2:'自适应多证据修复的信息流。上层由真实局部记录进入时间、波长与低秩背景三类证据，经持续性及参考一致性校准形成综合评分，再由连通跨度规则确定修复集合；下层按箭头依次进行秩4 Huber补全、稳定性检查、观测保真和局部非负支撑，最后由修复矩阵构造双线性连续函数。虚线独立表示顶码质量标记与条件性下界估计分支。案例数值来自6B03205在14.5 s、512.659 nm的正式记录和最终评分，图中流程不增加算法。',
6:'全部实际修复的时空分布及幅度。(a)逐点绘制1043个修复节点，横坐标为实际波长，纵向每个炉次带内由上至下表示相对冶炼时间0至1；颜色为有符号估计变化D−Y，圆点和菱形分别表示非回退与回退节点，右侧为逐炉修复数。(b)比较两子集绝对改变量的经验分布；(c)按各子集总数归一化展示相对时间分布。统计仅描述实际输出改动，不将节点作为独立炉次样本。',
7:'十炉连续光谱的紧凑曲面图。所有炉次采用相同正交投影、强度范围和共享色标，保留各炉实际波长范围及采集终点；各时间轴起点均为0。曲面显示网格包含全局峰值及已核实案例，关闭坐标区裁剪以保留完整前缘，并验证所有投影网格均位于整图范围内。独立大图及连续查询参数仍见既有附录和正式接口。'
}
def main():
 for r in pd.read_csv(O/'source_manifest.csv').itertuples():verify('source:'+r.path,hashlib.sha256((F.parents[1]/r.path).read_bytes()).hexdigest()==r.sha256)
 src=R/'问题一_中文论文_V3_研究深化与图表重构版.md';s=src.read_text();before=s
 s=s.replace('../figures/publication_v3/main/','../figures/publication_v4/main/')
 for n,c in caps.items():s=re.sub(rf'^图5-{n}　[^\n]+',f'图5-{n}　{c}',s,flags=re.M)
 s=re.sub(r'^为说明各项统计量如何作用于同一个观测，图5-2[^\n]+','图5-2将原始记录、证据构造、评分判定、低秩估计与连续输出连成一个有向流程。三类证据分别刻画时间偏离、局部谱形偏离和共同背景偏离，持续性与参考一致性调节其组合；修复集合只决定哪些节点进入估计，保留节点通过投影恢复原始细节。顶码观测沿独立分支保留记录量程信息，并形成条件性潜在估计。',s,flags=re.M)
 s=re.sub(r'^图5-6将逐炉修复规模与局部作用范围联系起来。[^\n]+','图5-6将全部修复节点放回各炉的波长—相对时间位置，展示修复在不同波段和阶段的聚集情况。图中的465个回退节点属于1043点内部子集，其绝对改变量中位数为16629.50 counts；其余578点对应221.36 counts。幅度分布与时间分布共同补充逐炉计数，使正式矩阵的改动能够按节点追溯。带内相对时间只用于跨炉绘图，连续函数仍使用实际秒坐标。',s,flags=re.M)
 s=re.sub(r'^图5-7汇总全部十炉的连续重构结果。[^\n]+','图5-7采用统一投影和强度尺度展示全部十炉结果，同时保留各炉采集终点。紧凑排列与共享单位减少重复标注，使曲面主体在整页中占据更充分的画幅；完整网格前缘、窄峰和代表案例均参与展示。各炉独立曲面见附图A1—A10，精确查询仍由式(5-23)定义。',s,flags=re.M)
 out=R/'问题一_中文论文_V4_机制与结果图优化版.md';out.write_text(s)
 verify('equations_unchanged',re.findall(r'\\tag\{[^}]+\}',before)==re.findall(r'\\tag\{[^}]+\}',s))
 verify('seven_sections',len(re.findall(r'^## 5\.',s,re.M))==7);verify('seven_figures',len(re.findall(r'^!\[',s,re.M))==7)
 for link in re.findall(r'\]\(([^)]+)\)',s):
  if not link.startswith(('http','mailto:')):verify('link:'+link,(R/link.split('#')[0]).exists())
 # Check each plotted node against the immutable interface, not merely another CSV.
 points=pd.read_csv(O/'repair_atlas_points.csv');verify('1043_unique_repair_ids',len(points)==1043 and not points.duplicated(['heat_id','row','column']).any())
 for h,d in points.groupby('heat_id'):
  with np.load(F/f'data/{h}_interface.npz') as z:
   i=d.row.to_numpy();j=d.column.to_numpy();verify('actual_mask:'+h,set(zip(i,j))==set(zip(*np.where(z['statistical_repair']))));verify('actual_difference:'+h,np.allclose(d.difference,z['reconstructed'][i,j]-z['raw'][i,j],rtol=0,atol=1e-8));verify('phase:'+h,np.allclose(d.phase,z['time'][i]/z['time'][-1]));verify('fallback:'+h,np.array_equal(d.numerical_fallback,z['numerical_fallback'][i,j]))
 c=json.loads((O/'checks.json').read_text());verify('frozen_partition',c['atlas']['fallback']==465 and c['atlas']['nonfallback']==578)
 for key in ['all_surface_maxima_retained','no_surface_axes_clipping','all_projected_meshes_within_figure']:verify(key,c[key])
 for p in ['F01_data_structure','F03_local_evidence','F04_continuous_reconstruction','F05_joint_validation']:
  for ext in ['png','svg','pdf']:verify('unchanged:'+p+ext,(O/f'main/{p}.{ext}').read_bytes()==(F/f'figures/publication_v3/main/{p}.{ext}').read_bytes())
 for row in pd.read_csv(O/'asset_manifest.csv').itertuples():
  p=O/row.path;verify('image_hash:'+row.path,hashlib.sha256(p.read_bytes()).hexdigest()==row.sha256)
  if row.format=='svg':verify('English:'+row.figure,not re.search(r'[\u3400-\u9fff]',p.read_text()))
  if row.format=='png':
   im=Image.open(p);verify('600dpi:'+row.figure,abs(im.info['dpi'][0]-600)<1)
   if row.figure.startswith(('F02','F06','F07')):
    im=im.convert('RGB');im.thumbnail((816,1220));page=Image.new('RGB',(992,1403),'white');page.paste(im,((992-im.width)//2,55));ImageDraw.Draw(page).text((88,20),row.figure+' | A4 page proof',fill='black');page.save(O/f'review/{row.figure}_page.png')
 # Confirm all previously published image bytes match committed files, despite reading them here.
 for p in (F/'figures/publication_v3').rglob('*'):
  if p.is_file() and p.suffix in ['.png','.svg','.pdf']:
   old=subprocess.run(['git','show',f'HEAD:{p.relative_to(F.parents[1])}'],capture_output=True,check=True).stdout;verify('v3_preserved:'+p.name,hashlib.sha256(old).hexdigest()==hashlib.sha256(p.read_bytes()).hexdigest())
 (O/'validation.json').write_text(json.dumps({'passed':len(checks),'failed':0,'model_or_data_changed':False,'original_v3_preserved':True,'changed_main_figures':['F02','F06','F07'],'other_main_figures_byte_identical':True},indent=2)+'\n')
 (O/'figure_captions.md').write_text('# 新版三图中文图注\n\n'+'\n\n'.join(f'图5-{n}　{c}' for n,c in caps.items())+'\n')
 print('PACKAGE_VERIFIED',len(checks))
if __name__=='__main__':main()
