"""Actual Times New Roman; centered panel titles; no scientific data changes."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
O=Path(__file__).resolve().parents[1]
FONTS=Path.home()/'.local/share/fonts/msttcore'
for p in FONTS.glob('*.TTF'):font_manager.fontManager.addfont(str(p))
actual=font_manager.findfont('Times New Roman',fallback_to_default=False)
plt.rcParams.update({'font.family':'Times New Roman','axes.titlelocation':'center','mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold','mathtext.cal':'Times New Roman','mathtext.sf':'Times New Roman','mathtext.tt':'Times New Roman','mathtext.fallback':'stix'})
RECORDS=json.loads((O/'typography_validation.json').read_text())['figures'] if (O/'typography_validation.json').exists() else {}
def digest(fig):
 h=hashlib.sha256()
 for ax in fig.axes:
  for a in ax.lines:
   for v in [a.get_xdata(),a.get_ydata()]:h.update(np.asarray(v).tobytes())
  for a in ax.collections:
   for v in [a.get_array(),a.get_offsets()]:
    if v is not None:h.update(np.asanyarray(v).tobytes())
  for a in ax.images:h.update(np.asanyarray(a.get_array()).tobytes())
 return h.hexdigest()
def apply(fig,name):
 before=digest(fig);titles=[]
 for ax in fig.axes:
  for old in [ax._left_title,ax._right_title]:
   if old.get_text():
    text=old.get_text();size=old.get_fontsize();weight=old.get_fontweight();old.set_text('');ax.set_title(text,loc='center',fontsize=size,fontweight=weight,pad=6)
  if ax.title.get_text():titles.append({'text':ax.title.get_text(),'x':ax.title.get_position()[0],'ha':ax.title.get_ha()})
 for t in fig.findobj(match=Text):t.set_fontfamily('Times New Roman')
 fig.canvas.draw();after=digest(fig);assert before==after
 assert all(x['x']==.5 and x['ha']=='center' for x in titles)
 RECORDS[name]={'font':'Times New Roman','centered_titles':titles,'plot_data_before':before,'plot_data_after':after,'plot_data_unchanged':before==after}
 (O/'typography_validation.json').write_text(json.dumps({'font_file':actual,'font_sha256':hashlib.sha256(Path(actual).read_bytes()).hexdigest(),'source':'https://downloads.sourceforge.net/corefonts/times32.exe','fonts_distributed_in_repository':False,'figures':RECORDS},indent=2)+'\n')
