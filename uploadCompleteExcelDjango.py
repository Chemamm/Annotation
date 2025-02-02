#!/usr/bin/env/python3

from Users.Arlong.Desktop.Doctorado.Tick_host import Annotation
from transcripts.models import Transcript

f = Annotation.CompleteExcel("/Users/Arlong/Desktop/Doctorado/Parasito/Transcriptome_assembly/MGSG_v2/MGSG_cd_hit_50aa_FPKM_v2_CompleteAnnotated_Unigenes.tsv")

for n in range(1,len(f.lines)):
    line = Annotation.CompleteExcel.Line(f.lines[n])
    if line.identity == "None":
        ident = 0
    else:
        ident = float(line.identity)
    if line.coverage == "None":
        cov = 0
    else:
        cov = int(line.coverage)
    if line.bitscore == "None":
        bit = 0
    else:
        bit = float(line.bitscore)
    if line.tm_domain == "None":
        tm_d = 0
    else:
        tm_d = int(line.tm_domain)
    if line.tm_domain_all == "None":
        tm_d_all = 0
    else:
        tm_d_all = int(line.tm_domain_all)
    foo = Transcript(cdsID = line.cdsID, transcriptID = line.transcriptID, transcript_seq = line.transcript_seq, cds_seq = line.cds_seq, pep_seq = line.seq, family = line.family, hitID = line.hitID, identity = ident, coverage = cov, bitscore = bit, description = line.description, db = line.db, GOterm = line.GOterm, GOterm_des = line.GOterm_des, KW = line.KW, KW_des = line.KW_des, ipro_id = line.ipro_id, ipro_des = line.ipro_des, ipro_goterm = line.ipro_goterm, ipro_goterm_des = line.ipro_goterm_des, merged = line.merged, merged_des = line.merged_des, sigp = line.sigp, tm_domain_all = tm_d_all, tm_domain = tm_d, classification = line.classification, eclass = line.eclass, SG_FC = float(line.SG_FC), MG_FC = float(line.MG_FC), total_FPKM = float(line.total_FPKM), DE_str = line.DE_str, DE_FC = line.DE_FC, DE_PPDE = line.DE_PPDE, SG_unfed_1_5mg_pre2_FPKM = float(line.row[32]), SG_unfed_1_7mg_pre4_FPKM = float(line.row[33]), SG_unfed_1_4mg_pre6_FPKM = float(line.row[34]), SG_unfed_1_6mg_pre8_FPKM = float(line.row[35]), SG_unfed_1_4mg_pre10_FPKM = float(line.row[36]), SG_1_12h_1_5mg_1_FPKM = float(line.row[37]), SG_1_12h_2mg_105_FPKM = float(line.row[38]), SG_1_12h_1_8mg_5_FPKM = float(line.row[39]), SG_1_12h_1_7mg_101_FPKM = float(line.row[40]), SG_1_12h_2mg_103_FPKM = float(line.row[41]), SG_1_24h_1_5mg_7_FPKM = float(line.row[42]), SG_1_24h_1_7mg_9_FPKM = float(line.row[43]), SG_1_24h_2_3mg_109_FPKM = float(line.row[44]), SG_1_24h_1_6mg_106_FPKM = float(line.row[45]), SG_1_24h_2_3mg_8_FPKM = float(line.row[46]), SG_1_48h_3_2mg_13_FPKM = float(line.row[47]), SG_1_48h_3mg_15_FPKM = float(line.row[48]), SG_1_48h_2_9mg_111_FPKM = float(line.row[49]), SG_1_48h_3_1mg_113_FPKM = float(line.row[50]), SG_1_48h_2_9mg_114_FPKM = float(line.row[51]), SG_1_72h_5_2mg_18_FPKM = float(line.row[52]), SG_1_72h_7_4mg_19_FPKM = float(line.row[53]), SG_1_72h_5_5mg_20_FPKM = float(line.row[54]), SG_1_72h_6_2mg_116_FPKM = float(line.row[55]), SG_1_72h_5_9mg_118_FPKM = float(line.row[56]), SG_1_96h_10_3mg_21_FPKM = float(line.row[57]), SG_1_96h_10_3mg_22_FPKM = float(line.row[58]), SG_1_96h_10_2mg_24_FPKM = float(line.row[59]), SG_1_96h_9_1mg_121_FPKM = float(line.row[60]), SG_1_96h_9_5mg_122_FPKM = float(line.row[61]), SG_2_12h_2_2mg_52_FPKM = float(line.row[62]), SG_2_12h_1_9mg_53_FPKM = float(line.row[63]), SG_2_12h_2mg_152_FPKM = float(line.row[64]), SG_2_12h_1_7mg_153_FPKM = float(line.row[65]), SG_2_12h_1_8mg_155_FPKM = float(line.row[66]), SG_2_24h_2_2mg_56_FPKM = float(line.row[67]), SG_2_24h_2_4mg_59_FPKM = float(line.row[68]), SG_2_24h_2_2mg_156_FPKM = float(line.row[69]), SG_2_24h_2_3mg_157_FPKM = float(line.row[70]), SG_2_24h_2_5mg_159_FPKM = float(line.row[71]), SG_2_48h_4_6mg_62_FPKM = float(line.row[72]), SG_2_48h_4_5mg_63_FPKM = float(line.row[73]), SG_2_48h_4_2mg_65_FPKM = float(line.row[74]), SG_2_48h_5mg_162_FPKM = float(line.row[75]), SG_2_48h_4_3mg_163_FPKM = float(line.row[76]), SG_2_72h_8_7mg_67_FPKM = float(line.row[77]), SG_2_72h_8_2mg_68_FPKM = float(line.row[78]), SG_2_72h_9_3mg_69_FPKM = float(line.row[79]), SG_2_72h_6_7mg_166_FPKM = float(line.row[80]), SG_2_72h_7_2mg_167_FPKM = float(line.row[81]), SG_2_96h_11_1mg_71_FPKM = float(line.row[82]), SG_2_96h_12_9mg_74_FPKM = float(line.row[83]), SG_2_96h_15_1mg_75_FPKM = float(line.row[84]), SG_2_96h_10_2mg_171_FPKM = float(line.row[85]), SG_2_96h_10_1mg_172_FPKM = float(line.row[86]), MG_unfed_1_7mg_pre4_FPKM = float(line.row[87]), MG_unfed_1_4mg_pre10_FPKM = float(line.row[88]), MG_unfed_1_6mg_pre8_FPKM = float(line.row[89]), MG_1_12h_1_8mg_4_FPKM = float(line.row[90]), MG_1_12h_1_8mg_5_FPKM = float(line.row[91]), MG_1_12h_1_8mg_104_FPKM = float(line.row[92]), MG_1_24h_1_7mg_9_FPKM = float(line.row[93]), MG_1_24h_1_6mg_10_FPKM = float(line.row[94]), MG_1_24h_1_9mg_107_FPKM = float(line.row[95]), MG_1_48h_3_2mg_13_FPKM = float(line.row[96]), MG_1_48h_3mg_15_FPKM = float(line.row[97]), MG_1_48h_2_9mg_114_FPKM = float(line.row[98]), MG_1_72h_5_5mg_20_FPKM = float(line.row[99]), MG_1_72h_6_2mg_116_FPKM = float(line.row[100]), MG_1_72h_5_9mg_118_FPKM = float(line.row[101]), MG_1_96h_10_3mg_22_FPKM = float(line.row[102]), MG_1_96h_10_2mg_24_FPKM = float(line.row[103]), MG_1_96h_9_5mg_122_FPKM = float(line.row[104]), MG_2_12h_1_8mg_155_FPKM = float(line.row[105]), MG_2_12h_2_2mg_54_FPKM = float(line.row[106]), MG_2_12h_2mg_152_FPKM = float(line.row[107]), MG_2_24h_2_2mg_56_FPKM = float(line.row[108]), MG_2_24h_2_4mg_59_FPKM = float(line.row[109]), MG_2_24h_2_3mg_157_FPKM = float(line.row[110]), MG_2_48h_4_6mg_62_FPKM = float(line.row[111]), MG_2_48h_4_5mg_63_FPKM = float(line.row[112]), MG_2_48h_4_3mg_163_FPKM = float(line.row[113]), MG_2_72h_8_7mg_67_FPKM = float(line.row[114]), MG_2_72h_8_2mg_68_FPKM = float(line.row[115]), MG_2_72h_7_2mg_167_FPKM = float(line.row[116]), MG_2_96h_11_1mg_71_FPKM = float(line.row[117]), MG_2_96h_12_9mg_74_FPKM = float(line.row[118]), MG_2_96h_10_2mg_171_FPKM = float(line.row[119]))
        foo.save()
