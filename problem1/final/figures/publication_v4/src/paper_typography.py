"""Review existing figures with external, editable English paper captions."""
from pathlib import Path
import json, hashlib, subprocess, xml.etree.ElementTree as ET
from PIL import Image,ImageDraw,ImageFont
O=Path(__file__).resolve().parents[1];ROOT=O.parents[3];R=O.parents[1]/'reports'
CAPTIONS={
'F01_data_structure':'Figure 5-1. Spectral structure and evidence for low-rank representation.',
'F02_model_mechanism':'Figure 5-2. Adaptive evidence fusion and robust low-rank spectral reconstruction.',
'F03_local_evidence':'Figure 5-3. Isolated anomaly repair and persistent spectral feature preservation.',
'F04_continuous_reconstruction':'Figure 5-4. Interpolation errors under time and wavelength downsampling.',
'F05_joint_validation':'Figure 5-5. Joint validation of anomaly detection and spectral repair.',
'F06_actual_repair_results':'Figure 5-6. Distribution and magnitude of the 1,043 repaired observations.',
'F07_ten_heat_surfaces':'Figure 5-7. Reconstructed continuous spectral surfaces for all ten heats.'}
def main():
 evidence=[];tmp=Path('/tmp/apmcm-paper-style');tmp.mkdir(exist_ok=True)
 for year,page,printed in [(2022,6,4),(2023,8,6),(2024,14,14),(2025,13,10)]:
  p=next((ROOT/'往年的优秀论文').rglob(f'{year}年B题优秀论文.pdf'));xml=tmp/f'{year}.xml'
  subprocess.run(['pdftohtml','-xml','-hidden','-i',str(p),str(xml)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
  doc=ET.parse(xml).getroot();fonts={x.attrib['id']:x.attrib for x in doc.iter('fontspec')};pg=next(x for x in doc.findall('page') if int(x.attrib['number'])==page)
  spans=[]
  for t in pg.findall('text'):
   text=''.join(t.itertext()).strip()
   if text.startswith(('Figure','Fig.')):
    f=fonts[t.attrib['font']];spans.append({'text':text,'font':f['family'],'size_pt_approx':round(float(f['size'])/1.5,2)})
  evidence.append({'year':year,'source':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'pdf_page':page,'printed_page':printed,'caption_spans':spans,'placement':'below figure, outside plot image','in_plot_font':'2022 vector mathematical labels: CMMI10/CMR7 and TeXGyreTermesX; other selected chart labels are not identifiable from extracted text; family not asserted'})
 (O/'paper_style_sources.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n')
 font=ImageFont.truetype('/usr/share/fonts/truetype/liberation2/LiberationSerif-Bold.ttf',18) # 10.8 pt at 120 dpi
 for name in ['F02_model_mechanism','F03_local_evidence','F04_continuous_reconstruction']:
  im=Image.open(O/f'main/{name}.png').convert('RGB');im.thumbnail((816,1120));page=Image.new('RGB',(992,1403),'white');y=100;page.paste(im,((992-im.width)//2,y));draw=ImageDraw.Draw(page);words=CAPTIONS[name].split();lines=[];line=''
  for word in words:
   trial=(line+' '+word).strip()
   if draw.textlength(trial,font=font)>816:lines.append(line);line=word
   else:line=trial
  lines.append(line)
  for k,line in enumerate(lines):draw.text((496,y+im.height+18+k*24),line,font=font,anchor='mt',fill='black')
  page.save(O/f'review/{name}_caption_page.png',dpi=(120,120))
 (O/'english_captions.md').write_text('# English figure captions for final typesetting\n\nCaptions are editable manuscript text placed below the figure, outside PNG/SVG/PDF assets.\n\n'+'\n\n'.join(CAPTIONS.values())+'\n')
if __name__=='__main__':main()
