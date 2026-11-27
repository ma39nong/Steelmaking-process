# 光谱物理参考与波长口径

核读日期：2026-11-27。以下来源用于解释候选结构，附件实际波长轴和统计事实另由计算结果给出。

|来源|核读内容|用于本题的范围|
|---|---|---|
|[NIST ASD](https://www.nist.gov/pml/atomic-spectra-database)、[谱线查询说明](https://physics.nist.gov/PhysRefData/ASD/lines_form.html)|数据库及空气/真空输出选项已核读；动态Na查询入口本次未成功返回|采用下面NIST官方强线表可核实的参考值；不伪称动态查询取得新结果|
|[NIST Sodium strong lines](https://physics.nist.gov/PhysRefData/Handbook/Tables/sodiumtable2.htm)|Air wavelength栏的Na I 5889.950、5895.924 Å|参考空气波长588.9950、589.5924 nm；约0.5974 nm间距|
|[NIST Potassium strong lines](https://physics.nist.gov/PhysRefData/Handbook/Tables/potassiumtable2.htm)|Air wavelength栏的K I 7664.8991、7698.9645 Å|参考空气波长766.48991、769.89645 nm；约3.40654 nm间距|
|Brown L R, Plymate C. Experimental Line Parameters of the Oxygen A Band at 760 nm. Journal of Molecular Spectroscopy, 2000, 199(2):166–179. DOI10.1006/jmsp.1999.8012，[PubMed](https://pubmed.ncbi.nlm.nih.gov/10637102/)|书目与摘要已读，未取得全文|760 nm附近的吸收带背景；不作为附件已识别O₂吸收的证据|
|Ocean Optics USB4000 manual，Document211-00000-000-02-201604，[官方PDF](https://www.oceanoptics.com/wp-content/uploads/2024/12/USB4000-Operating-Instructions.pdf)|已有合法全文`references/fulltext/R09.pdf`，PDF第10页列出16位ADC与0.1–10 nm FWHM范围，明确取决于光栅和狭缝|附件对应的实际仪器配置未知，不能以通道间隔代替光学分辨率|
|Zhao等，Scientific Reports13:14409，2023，DOI10.1038/s41598-023-41751-9，[PMC入口](https://pmc.ncbi.nlm.nih.gov/articles/PMC10475112/)|本次网页受访问验证限制；沿用既有已核读OA全文`references/fulltext/R04_0.pdf`，pp.2–4|转炉采集及预处理背景，与附件网格和仪器配置分别表达|

上述原子参考是空气波长。附件表头未注明空气/真空约定，因而`physical_feature_candidates.csv`中的差值仅为实际通道与指定参考坐标的距离，不解释为仪器标定误差。小于通道间距的接近程度可支持位置一致性，但不能单独确定物种归属。Na双线是否可分辨还取决于未知的光学FWHM与谱形；K双线及760 nm附近结构同样需结合邻域、时间持续性和亮度关联。

6C02108在766.413 nm的保留通道与K I空气参考766.48991 nm相距−0.07691 nm。该点有持续谱形证据和位置一致性，保留为物理候选特征；统计事实与元素鉴定分别记录。383.293、835.887 nm作为持续窄结构诊断目标，不直接命名为故障或化学组分。
