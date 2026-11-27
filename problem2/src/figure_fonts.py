"""Shared strict Times New Roman configuration for new Q2 audit figures."""
from pathlib import Path
import hashlib
import subprocess
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager as fm


def configure():
    # Fontconfig can know recently installed fonts absent from Matplotlib's cache.
    for face in ['Regular','Bold','Italic','Bold Italic']:
        path=subprocess.check_output(['fc-match','-f','%{file}','Times New Roman:style='+face],text=True).strip()
        if fm.FontProperties(fname=path).get_name()!='Times New Roman':
            raise RuntimeError('real Times New Roman is required')
        fm.fontManager.addfont(path)
    records=[]
    for style,weight in [('normal','normal'),('normal','bold'),('italic','normal'),('italic','bold')]:
        prop=fm.FontProperties(family='Times New Roman',style=style,weight=weight)
        path=Path(fm.findfont(prop,fallback_to_default=False))
        family=fm.FontProperties(fname=path).get_name()
        if family!='Times New Roman': raise RuntimeError('Times New Roman unavailable; refusing substitution')
        records.append(dict(style=style,weight=weight,family=family,file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    matplotlib.rcParams.update({'font.family':'Times New Roman','font.size':9,
        'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
        'legend.fontsize':8,'svg.hashsalt':'problem2-delay-audit','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'mathtext.fontset':'custom','mathtext.rm':'Times New Roman',
        'mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold',
        'mathtext.fallback':'stix','axes.unicode_minus':False})
    return dict(fonts=records,special_math_fallback='STIX explicitly configured; audit figures use plain Latin text and digits',
                intended_width_inches=7.2,svg_text='editable; requires installed Times New Roman')
