#!/usr/bin/env/python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pylab import *
import os
from os.path import exists
import math
from scipy.stats import spearmanr, variation
import math
import Plotting
import random
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import MDS
import scipy.cluster.hierarchy as shc
from matplotlib_venn import venn2
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.obo_parser import GODag
from matplotlib.lines import Line2D
import scipy



class CompleteExcel:
    def __init__(self, file):
        with open(file) as ex:
            self.lines = list(ex)
        self.IDs, self.transcriptIDs, self.IDs_GO, self.IDs_KW, self.IDs_IP, self.IDs_family, self.IDs_merged = self.getIDs()
        self.df = pd.read_csv(file, sep="\t", index_col=0)
        self.DE_dict = self.getDE()
        self.annotated = union(self.IDs_merged, self.IDs_IP)
        self.GO_dict, self.GO_description_dict = self.getGO()
        self.secreted = self.getSecreted()
        self.uniprot_dict = self.get_uniprot_dict()


    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.cdsID = self.row[0]
            self.transcriptID = self.cdsID.split(".")[0]
            self.transcript_seq = self.row[1]
            self.cds_seq = self.row[2]
            self.seq = self.row[3]
            self.family = self.row[4]
            self.hitID = self.row[5]
            self.identity = self.row[6]
            self.coverage = self.row[7]
            self.bitscore = self.row[8]
            self.description = self.row[9]
            self.db = self.row[10]
            self.GOterm = self.row[11]
            self.GOterm_des = self.row[12]
            self.KW = self.row[13]
            self.KW_des = self.row[14]
            self.ipro_id = self.row[15]
            self.ipro_des = self.row[16]
            self.ipro_goterm = self.row[17]
            self.ipro_goterm_des = self.row[18]
            self.merged = self.row[19]
            self.merged_des = self.row[20]
            self.sigp = self.row[21]
            self.tm_domain_all = self.row[22]
            self.tm_domain = self.row[23]
            self.classification = self.row[24]
            self.eclass = self.row[25]
            self.SG_FC = self.row[26]
            self.MG_FC = self.row[27]
            self.total_FPKM = self.row[28]
            self.DE_str = self.row[29]
            self.DE_FC = self.row[30]
            self.DE_PPDE = self.row[31]
            self.evalues = self.row[32::]

    def getIDs(self):
        IDs = []
        transcriptIDs = []
        IDs_GO = []
        IDs_KW = []
        IDs_IP = []
        IDs_family = []
        IDs_merged = []
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            IDs.append(line.cdsID)
            transcriptIDs.append(line.transcriptID)
            if line.GOterm != "None":
                IDs_GO.append(line.transcriptID)
            if line.KW != "None":
                IDs_KW.append(line.transcriptID)
            if line.ipro_id != "None":
                IDs_IP.append(line.transcriptID)
            if line.family != "None":
                IDs_family.append(line.transcriptID)
            if line.merged != "None":
                IDs_merged.append(line.transcriptID)
        return IDs, transcriptIDs, IDs_GO, IDs_KW, IDs_IP, IDs_family, IDs_merged

    def getGO(self):
        GO_dict = dict()
        GO_description_dict = dict()
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.merged != "None":
                GO_dict[line.transcriptID] = line.merged.split("; ")
                GO_description_dict[line.transcriptID] = line.merged_des.split("; ")
        return GO_dict, GO_description_dict

    def getGOBackground(self, output):
        with open(output, "wt") as out:
            for ID in self.GO_dict:
                out.write("%s\t%s\n" %(ID, ",".join(self.GO_dict[ID])))

    def get_uniprot_dict(self):
        uniprot_dict = dict()
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            uniprot_dict[line.transcriptID] = line.hitID
        return uniprot_dict



    def getSecreted(self):
        secreted = []
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.classification.startswith("Putatively Secreted"):
                secreted.append(line.transcriptID)
        return secreted

    def get_annot_file(self,output):
        GO_dict = {}
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.merged != "None":
                GO_dict[line.transcriptID] = line.merged.split("; ")
        with open(output,"wt") as out:
            for key in GO_dict:
                out.write("%s\t%s\n" %(key, ",".join(GO_dict[key])))


    def splitSecreted(self, output_label):
        with open("%s_secreted.fa" %output_label, "wt") as secreted_out:
            with open("%s_non-secreted.fa" %output_label, "wt") as nosecreted_out:
                for n in range(1,len(self.lines)):
                    line = CompleteExcel.Line(self.lines[n])
                    if line.classification == "Putatively Secreted, Annotated" or line.classification == "Putatively Secreted, not Annotated":
                        secreted_out.write(">%s\n%s\n" %(line.cdsID,line.seq))
                    else:
                        nosecreted_out.write(">%s\n%s\n" %(line.cdsID,line.seq))
    def getDE(self):
        DE_dict = {}
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            DE_dict[line.cdsID] = line.DE_str.split("; ")
        return DE_dict

    def getPeptidefa(self, list, output):
        with open(output, "wt") as out:
            for n in range(1,len(self.lines)):
                line = CompleteExcel.Line(self.lines[n])
                if line.transcriptID in list:
                    out.write(">%s\n%s\n" %(line.transcriptID, line.seq))

    def orderbyweight(self):
        head_evalues = self.lines[0].replace("\n","").split("\t")[32::]
        MG_1_dict = {}
        MG_2_dict = {}
        SG_1_dict = {}
        SG_2_dict = {}
        for sample in head_evalues:
            tissue = sample.split("_")[0]
            exposure = sample.split("_")[1]
            weight = float(sample.split("mg")[0].split("_")[-1].replace(",","."))
            if tissue == "MG" and exposure == "1":
                MG_1_dict[sample] = weight
            elif tissue == "MG" and exposure == "2":
                MG_2_dict[sample] = weight
            elif tissue == "SG" and exposure == "1":
                SG_1_dict[sample] = weight
            elif tissue == "SG" and exposure == "2":
                SG_2_dict[sample] = weight
        #Ordering
        MG_1_dict = {k: v for k, v in sorted(MG_1_dict.items(), key=lambda item: item[1])}
        MG_2_dict = {k: v for k, v in sorted(MG_2_dict.items(), key=lambda item: item[1])}
        SG_1_dict = {k: v for k, v in sorted(SG_1_dict.items(), key=lambda item: item[1])}
        SG_2_dict = {k: v for k, v in sorted(SG_2_dict.items(), key=lambda item: item[1])}
        MG_1_order = list(MG_1_dict.keys())
        MG_2_order = list(MG_2_dict.keys())
        SG_1_order = list(SG_1_dict.keys())
        SG_2_order = list(SG_2_dict.keys())
        return MG_1_order, MG_2_order, SG_1_order, SG_2_order

    def getstats(self, output):
        GO = 0
        notGO = 0
        GO_SA = 0
        GO_SNA = 0
        GO_NSA = 0
        GO_NSNA = 0
        notGO_SA = 0
        notGO_SNA = 0
        notGO_NSA = 0
        notGO_NSNA = 0
        GO_interpro = 0
        notGO_interpro = 0
        GO_SA_SG = 0
        GO_SA_MG = 0
        GO_SA_N = 0
        GO_SNA_SG = 0
        GO_SNA_MG = 0
        GO_SNA_N = 0
        GO_NSA_SG = 0
        GO_NSA_MG = 0
        GO_NSA_N = 0
        GO_NSNA_SG = 0
        GO_NSNA_MG = 0
        GO_NSNA_N = 0
        notGO_SA_SG = 0
        notGO_SA_MG = 0
        notGO_SA_N = 0
        notGO_SNA_SG = 0
        notGO_SNA_MG = 0
        notGO_SNA_N = 0
        notGO_NSA_SG = 0
        notGO_NSA_MG = 0
        notGO_NSA_N = 0
        notGO_NSNA_SG = 0
        notGO_NSNA_MG = 0
        notGO_NSNA_N = 0
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.GOterm != "None":
                GO += 1
                if line.classification == "Putatively Secreted, Annotated":
                    GO_SA += 1
                    if line.eclass == "SG specific":
                        GO_SA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_SA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    GO_SNA += 1
                    if line.eclass == "SG specific":
                        GO_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    GO_NSA += 1
                    if line.eclass == "SG specific":
                        GO_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    GO_NSNA += 1
                    if line.eclass == "SG specific":
                        GO_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_NSNA_N += 1
                if line.ipro_id != "None":
                    GO_interpro += 1

            else:
                notGO += 1
                if line.classification == "Putatively Secreted, Annotated":
                    notGO_SA += 1
                    if line.eclass == "SG specific":
                        notGO_SA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_SA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    notGO_SNA += 1
                    if line.eclass == "SG specific":
                        notGO_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    notGO_NSA += 1
                    if line.eclass == "SG specific":
                        notGO_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    notGO_NSNA += 1
                    if line.eclass == "SG specific":
                        notGO_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_NSNA_N += 1
                if line.ipro_id != "None":
                    notGO_interpro += 1

        with open("%s_uniprotGO.stats" % output, "wt") as out:
            out.write("Total of CDS with GO annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (GO, GO_interpro, GO_SA, GO_SA_SG, GO_SA_MG, GO_SA_N, GO_SNA, GO_SNA_SG, GO_SNA_MG, GO_SNA_N, GO_NSA, GO_NSA_SG, GO_NSA_MG, GO_NSA_N, GO_NSNA, GO_NSNA_SG, GO_NSNA_MG, GO_NSNA_N))
            out.write("Total of CDS without GO annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (notGO, notGO_interpro, notGO_SA, notGO_SA_SG, notGO_SA_MG, notGO_SA_N, notGO_SNA, notGO_SNA_SG, notGO_SNA_MG, notGO_SNA_N, notGO_NSA, notGO_NSA_SG, notGO_NSA_MG, notGO_NSA_N, notGO_NSNA, notGO_NSNA_SG, notGO_NSNA_MG, notGO_NSNA_N))

        #Interpro
        ipro = 0
        notipro = 0
        ipro_SA = 0
        ipro_SNA = 0
        ipro_NSA = 0
        ipro_NSNA = 0
        notipro_SA = 0
        notipro_SNA = 0
        notipro_NSA = 0
        notipro_NSNA = 0
        ipro_GO = 0
        notipro_GO = 0
        ipro_SA_SG = 0
        ipro_SA_MG = 0
        ipro_SA_N = 0
        ipro_SNA_SG = 0
        ipro_SNA_MG = 0
        ipro_SNA_N = 0
        ipro_NSA_SG = 0
        ipro_NSA_MG = 0
        ipro_NSA_N = 0
        ipro_NSNA_SG = 0
        ipro_NSNA_MG = 0
        ipro_NSNA_N = 0
        notipro_SA_SG = 0
        notipro_SA_MG = 0
        notipro_SA_N = 0
        notipro_SNA_SG = 0
        notipro_SNA_MG = 0
        notipro_SNA_N = 0
        notipro_NSA_SG = 0
        notipro_NSA_MG = 0
        notipro_NSA_N = 0
        notipro_NSNA_SG = 0
        notipro_NSNA_MG = 0
        notipro_NSNA_N = 0
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.ipro_id != "None":
                ipro += 1
                if line.classification == "Putatively Secreted, Annotated":
                    ipro_SA += 1
                    if line.eclass == "SG specific":
                        ipro_SA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_SA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    ipro_SNA += 1
                    if line.eclass == "SG specific":
                        ipro_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    ipro_NSA += 1
                    if line.eclass == "SG specific":
                        ipro_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    ipro_NSNA += 1
                    if line.eclass == "SG specific":
                        ipro_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_NSNA_N += 1
                if line.GOterm != "None":
                    ipro_GO += 1

            else:
                notipro += 1
                if line.classification == "Putatively Secreted, Annotated":
                    notipro_SA += 1
                    if line.eclass == "SG specific":
                        notipro_SA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_SA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    notipro_SNA += 1
                    if line.eclass == "SG specific":
                        notipro_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    notipro_NSA += 1
                    if line.eclass == "SG specific":
                        notipro_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    notipro_NSNA += 1
                    if line.eclass == "SG specific":
                        notipro_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_NSNA_N += 1
                if line.GOterm != "None":
                    notipro_GO += 1

        with open("%s_interpro.stats" % output, "wt") as out:
            out.write("Total of CDS with ipro annotation: %s. (%s of them have GO terms)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (ipro, ipro_GO, ipro_SA, ipro_SA_SG, ipro_SA_MG, ipro_SA_N, ipro_SNA, ipro_SNA_SG, ipro_SNA_MG, ipro_SNA_N, ipro_NSA, ipro_NSA_SG, ipro_NSA_MG, ipro_NSA_N, ipro_NSNA, ipro_NSNA_SG, ipro_NSNA_MG, ipro_NSNA_N))
            out.write("Total of CDS without ipro annotation: %s. (%s of them have GO terms)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (notipro, notipro_GO, notipro_SA, notipro_SA_SG, notipro_SA_MG, notipro_SA_N, notipro_SNA, notipro_SNA_SG, notipro_SNA_MG, notipro_SNA_N, notipro_NSA, notipro_NSA_SG, notipro_NSA_MG, notipro_NSA_N, notipro_NSNA, notipro_NSNA_SG, notipro_NSNA_MG, notipro_NSNA_N))
    #Merged GO
        merged = 0
        notmerged = 0
        merged_SA = 0
        merged_SNA = 0
        merged_NSA = 0
        merged_NSNA = 0
        notmerged_SA = 0
        notmerged_SNA = 0
        notmerged_NSA = 0
        notmerged_NSNA = 0
        merged_interpro = 0
        notmerged_interpro = 0
        merged_SA_SG = 0
        merged_SA_MG = 0
        merged_SA_N = 0
        merged_SNA_SG = 0
        merged_SNA_MG = 0
        merged_SNA_N = 0
        merged_NSA_SG = 0
        merged_NSA_MG = 0
        merged_NSA_N = 0
        merged_NSNA_SG = 0
        merged_NSNA_MG = 0
        merged_NSNA_N = 0
        notmerged_SA_SG = 0
        notmerged_SA_MG = 0
        notmerged_SA_N = 0
        notmerged_SNA_SG = 0
        notmerged_SNA_MG = 0
        notmerged_SNA_N = 0
        notmerged_NSA_SG = 0
        notmerged_NSA_MG = 0
        notmerged_NSA_N = 0
        notmerged_NSNA_SG = 0
        notmerged_NSNA_MG = 0
        notmerged_NSNA_N = 0
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.merged != "None":
                merged += 1
                if line.classification == "Putatively Secreted, Annotated":
                    merged_SA += 1
                    if line.eclass == "SG specific":
                        merged_SA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_SA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    merged_SNA += 1
                    if line.eclass == "SG specific":
                        merged_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    merged_NSA += 1
                    if line.eclass == "SG specific":
                        merged_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    merged_NSNA += 1
                    if line.eclass == "SG specific":
                        merged_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_NSNA_N += 1
                if line.ipro_id != "None":
                    merged_interpro += 1

            else:
                notmerged += 1
                if line.classification == "Putatively Secreted, Annotated":
                    notmerged_SA += 1
                    if line.eclass == "SG specific":
                        notmerged_SA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_SA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    notmerged_SNA += 1
                    if line.eclass == "SG specific":
                        notmerged_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    notmerged_NSA += 1
                    if line.eclass == "SG specific":
                        notmerged_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    notmerged_NSNA += 1
                    if line.eclass == "SG specific":
                        notmerged_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_NSNA_N += 1
                if line.ipro_id != "None":
                    notmerged_interpro += 1

        with open("%s_merged.stats" % output, "wt") as out:
            out.write("Total of CDS with merged annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (merged, merged_interpro, merged_SA, merged_SA_SG, merged_SA_MG, merged_SA_N, merged_SNA, merged_SNA_SG, merged_SNA_MG, merged_SNA_N, merged_NSA, merged_NSA_SG, merged_NSA_MG, merged_NSA_N, merged_NSNA, merged_NSNA_SG, merged_NSNA_MG, merged_NSNA_N))
            out.write("Total of CDS without merged annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (notmerged, notmerged_interpro, notmerged_SA, notmerged_SA_SG, notmerged_SA_MG, notmerged_SA_N, notmerged_SNA, notmerged_SNA_SG, notmerged_SNA_MG, notmerged_SNA_N, notmerged_NSA, notmerged_NSA_SG, notmerged_NSA_MG, notmerged_NSA_N, notmerged_NSNA, notmerged_NSNA_SG, notmerged_NSNA_MG, notmerged_NSNA_N))

    def length_distribution(self, output):
        GOseq_dict = {}
        notGOseq_dict = {}
        iproseq_dict = {}
        notiproseq_dict = {}
        mergedseq_dict = {}
        notmergedseq_dict = {}
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.GOterm != "None":
                GOseq_dict[line.cdsID] = line.seq
            else:
                notGOseq_dict[line.cdsID] = line.seq
            if line.ipro_id != "None":
                iproseq_dict[line.cdsID] = line.seq
            else:
                notiproseq_dict[line.cdsID] = line.seq
            if line.merged != "None":
                mergedseq_dict[line.cdsID] = line.seq
            else:
                notmergedseq_dict[line.cdsID] = line.seq

        GOlength_dict = getlength(GOseq_dict)
        notGOlength_dict = getlength(notGOseq_dict)
        iprolength_dict = getlength(iproseq_dict)
        notiprolength_dict = getlength(notiproseq_dict)
        mergedlength_dict = getlength(mergedseq_dict)
        notmergedlength_dict = getlength(notmergedseq_dict)

        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(list(GOlength_dict.values()))
        ax.set(xlabel='With Uniprot/B2GO GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        subplot(1,2,2)
        ax = sns.histplot(list(notGOlength_dict.values()))
        ax.set(xlabel='Without Uniprot/B2GO GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        fig = ax.get_figure()
        fig.savefig("%s_uniprotGO_length_distribution.pdf" %output)
        plt.close()

        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(list(iprolength_dict.values()))
        ax.set(xlabel='With Interpro Annotation')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        subplot(1,2,2)
        ax = sns.histplot(list(notiprolength_dict.values()))
        ax.set(xlabel='Without Interpro Annotation')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        fig = ax.get_figure()
        fig.savefig("%s_interpro_length_distribution.pdf" %output)
        plt.close()

        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(list(mergedlength_dict.values()))
        ax.set(xlabel='With GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        subplot(1,2,2)
        ax = sns.histplot(list(notmergedlength_dict.values()))
        ax.set(xlabel='Without GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        fig = ax.get_figure()
        fig.savefig("%s_mergedGO_length_distribution.pdf" %output)
        plt.close()

    def expression_distribution(self, output):
        self.evalues_secreted, self.evalues_nonsecreted = get_expression_dicts(self.lines)
        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(data=self.evalues_secreted, log_scale=True)
        ax.set(xlabel='Secreted FPKM')
        subplot(1,2,2)
        ax = sns.histplot(data=self.evalues_nonsecreted, log_scale=True)
        ax.set(xlabel='Non Secreted FPKM')
        fig = ax.get_figure()
        fig.savefig(output)
        plt.close()

    def getBackgrounds(self,output_label):
        GO_dict = {}
        KW_dict = {}
        IP_dict = {}
        family_dict = {}
        classification_dict = {}
        tissue_dict = {}
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            #GO
            if line.merged != "None":
                GO_dict[line.transcriptID] = line.merged.replace("; ",",")
            #KW
            if line.KW != "None":
                KW_dict[line.transcriptID] = line.KW.replace("; ",",")
            #IP
            if line.ipro_id != "None":
                IP_dict[line.transcriptID] = line.ipro_id.replace("; ",",")
            #family
            if line.family != "None":
                family_dict[line.transcriptID] = line.family.replace("; ",",")
            #classification
            classification_dict[line.transcriptID] = line.classification.replace("; ",",")
            #classification
            tissue_dict[line.transcriptID] = line.eclass.replace("; ",",")

        with open("%s_GO_all.txt" %output_label,"wt") as out:
            for key in GO_dict:
                out.write("%s\t%s\n" %(key, GO_dict[key]))
        with open("%s_KW_all.txt" %output_label,"wt") as out:
            for key in KW_dict:
                out.write("%s\t%s\n" %(key, KW_dict[key]))
        with open("%s_IP_all.txt" %output_label,"wt") as out:
            for key in IP_dict:
                out.write("%s\t%s\n" %(key, IP_dict[key]))
        with open("%s_family_all.txt" %output_label,"wt") as out:
            for key in family_dict:
                out.write("%s\t%s\n" %(key, family_dict[key]))
        with open("%s_classification_all.txt" %output_label,"wt") as out:
            for key in classification_dict:
                out.write("%s\t%s\n" %(key, classification_dict[key]))
        with open("%s_eclass_all.txt" %output_label,"wt") as out:
            for key in tissue_dict:
                out.write("%s\t%s\n" %(key, tissue_dict[key]))

    def getannotationdicts(self):
        go_dict = {}
        kw_dict = {}
        ip_dict = {}
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            go = line.GOterm.split("; ")
            go_des = line.GOterm_des.split("; ")
            ip_go = line.ipro_goterm.split("; ")
            ip_go_des = line.ipro_goterm_des.split("; ")
            ip = line.ipro_id.split("; ")
            ip_des = line.ipro_des.split(";")
            kw = line.KW.split("; ")
            kw_des = line.KW_des.split(";")
            for n in range(0,len(go)):
                if len(go) == len(go_des):
                    go_dict[go[n]] = go_des[n]
                else:
                    print(go)
                    print(go_des)
            for n in range(0, len(ip_go)):
                go_dict[ip_go[n]] = ip_go_des[n]
            for n in range(0, len(kw)):
                kw_dict[kw[n]] = kw_des[n]
            for n in range(0,len(ip)):
                ip_dict[ip[n]] = ip_des[n]
        return go_dict, kw_dict, ip_dict

    def getUniprotIDs(self, gene_list, output):
        id_dict = dict()
        n = 0
        k = 0
        for gene in gene_list:
            if "." in gene:
                gene = gene.split(".")[0]
            if gene in self.uniprot_dict:
                n += 1
                id_dict[gene] = self.uniprot_dict[gene]
            else:
                k += 1
                id_dict[gene] = "ncRNA"

        with open(output, "wt") as out:
            out.write("gene\tUniprot ID\n")
            for gene in id_dict:
                out.write("%s\t%s\n" %(gene, id_dict[gene]))
        print("cRNA: %i\nncRNA: %i" %(n,k))

class TranscriptList:
    def __init__(self, transcript_list, sep=False):
        if sep:
            with open(transcript_list) as tl:
                if sep != "\n":
                    tl = tl.read().replace("\n","")
                    self.transcripts = tl.split(sep)
                else:
                    self.transcripts = list(tl)
        else:
            self.transcripts = transcript_list

    def depurate(self):
        self.transcripts = [x.split(".")[0] for x in self.transcripts]

    def annotationAnalysis(self, cexcel, output_label):
        cexcel = CompleteExcel(cexcel)
        annotated = []
        unannotated = []
        secreted = []
        non_secreted = []
        GO_descriptions = []
        with open("%s_GO.txt" %output_label, "wt") as out:
            for transcript in self.transcripts:
                if transcript in cexcel.GO_dict:
                    out.write("%s\t%s\n" %(transcript, ",".join(cexcel.GO_dict[transcript])))
                    annotated.append(transcript)
                    GO_descriptions.extend(cexcel.GO_description_dict[transcript])
                else:
                    unannotated.append(transcript)
                if transcript in cexcel.secreted:
                    secreted.append(transcript)
                else:
                    non_secreted.append(transcript)
        GO_descriptions_distribution = freq(GO_descriptions)
        GO_descriptions_distribution = {k: v for k, v in sorted(GO_descriptions_distribution.items(),
                                                                key=lambda item: item[1], reverse=True)}
        with open("%s_stats.txt" %output_label, "wt") as out:
            out.write("Annotated\t%i\nNot annotated\t%i\nSecreted\t%i\nNot secreted\t%i\nTotal GO\t"
                      "%i\n" %(len(annotated), len(unannotated), len(secreted), len(non_secreted), len(GO_descriptions)))
        with open("%s_GOdistribution.txt" %output_label, "wt") as out:
            for GO in GO_descriptions_distribution:
                out.write("%s\t%i\n" %(GO, GO_descriptions_distribution[GO]))

        df_dict = dict()
        df_dict["GO terms"] = GO_descriptions
        df = pd.DataFrame.from_dict(df_dict)
        fig = plt.figure()
        sns.set_theme()
        sns.countplot(x="GO terms", order=GO_descriptions_distribution.keys(),  data=df)
        plt.xticks(rotation=40)
        fig.savefig("%s_GOdistribution.png" %output_label)
        plt.close()

        ### 20 most represented
        go_20_list = []
        for n in range(0,20):
            GO_description = list(GO_descriptions_distribution.keys())[n]
            go_20_list.extend([GO_description] * GO_descriptions_distribution[GO_description])

        df_dict = dict()
        df_dict["GO terms"] = go_20_list
        df = pd.DataFrame.from_dict(df_dict)
        fig = plt.figure(figsize=(10, 5))
        sns.set_theme()
        sns.countplot(y="GO terms", data=df, color="blue")
        fig.tight_layout()
        fig.savefig("%s_GOdistribution_20.png" % output_label)
        plt.close()

class Ematrix:
    def __init__(self, file):
        self.name = file
        with open(file) as ex:
            self.lines = list(ex)
        self.df = pd.read_csv(file, sep="\t", index_col=0)
        self.dict = self.matrix2dict()
        self.header = self.lines[0]
        # Rename transcriptome Ematrix
        rename_dict = dict()
        rename_dict["pre_2_SG_S235"] = "SG_unfed_1,5mg_pre2"
        rename_dict["pre_4_SG_S267"] = "SG_unfed_1,7mg_pre4"
        rename_dict["pre_6_SG_S236"] = "SG_unfed_1,4mg_pre6"
        rename_dict["pre_8_SG_S274"] = "SG_unfed_1,6mg_pre8"
        rename_dict["pre_10_SG_S237"] = "SG_unfed_1,4mg_pre10"
        rename_dict["1_12h_SG_S268"] = "SG_1_12h_1,5mg_1"
        rename_dict["105_12h_SG_S1"] = "SG_1_12h_2mg_105"
        rename_dict["5_12h_SG_S224"] = "SG_1_12h_1,8mg_5"
        rename_dict["101_12h_SG_S225"] = "SG_1_12h_1,7mg_101"
        rename_dict["103_12h_SG_S275"] = "SG_1_12h_2mg_103"
        rename_dict["7_24h_SG_S270"] = "SG_1_24h_1,5mg_7"
        rename_dict["9_24h_SG_S226"] = "SG_1_24h_1,7mg_9"
        rename_dict["109_24h_SG_S312"] = "SG_1_24h_2,3mg_109"
        rename_dict["106_24h_SG_S227"] = "SG_1_24h_1,6mg_106"
        rename_dict["8_24h_SG_S269"] = "SG_1_24h_2,3mg_8"
        rename_dict["13_48h_SG_S255"] = "SG_1_48h_3,2mg_13"
        rename_dict["15_48h_SG_S256"] = "SG_1_48h_3mg_15"
        rename_dict["111_48h_SG_S228"] = "SG_1_48h_2,9mg_111"
        rename_dict["113_48h_SG_S229"] = "SG_1_48h_3,1mg_113"
        rename_dict["114_48h_SG_S271"] = "SG_1_48h_2,9mg_114"
        rename_dict["18_72h_SG_S257"] = "SG_1_72h_5,2mg_18"
        rename_dict["19_72h_SG_S258"] = "SG_1_72h_7,4mg_19"
        rename_dict["20_72h_SG_S259"] = "SG_1_72h_5,5mg_20"
        rename_dict["116_72h_SG_S307"] = "SG_1_72h_6,2mg_116"
        rename_dict["118_72h_SG_S230"] = "SG_1_72h_5,9mg_118"
        rename_dict["21_96h_SG_S239"] = "SG_1_96h_10,3mg_21"
        rename_dict["22_96h_SG_S240"] = "SG_1_96h_10,3mg_22"
        rename_dict["24_96h_SG_S260"] = "SG_1_96h_10,2mg_24"
        rename_dict["121_96h_SG_S231"] = "SG_1_96h_9,1mg_121"
        rename_dict["122_96h_SG_S261"] = "SG_1_96h_9,5mg_122"
        rename_dict["52_12h_SG_S280"] = "SG_2_12h_2,2mg_52"
        rename_dict["53_12h_SG_S232"] = "SG_2_12h_1,9mg_53"
        rename_dict["152_12h_SG_S233"] = "SG_2_12h_2mg_152"
        rename_dict["153_12h_SG_S234"] = "SG_2_12h_1,7mg_153"
        rename_dict["155_12h_SG_S272"] = "SG_2_12h_1,8mg_155"
        rename_dict["56_1d_SG_S262"] = "SG_2_24h_2,2mg_56"
        rename_dict["59_1d_SG_S263"] = "SG_2_24h_2,4mg_59"
        rename_dict["156_1d_SG_S241"] = "SG_2_24h_2,2mg_156"
        rename_dict["157_1d_SG_S273"] = "SG_2_24h_2,3mg_157"
        rename_dict["159_1d_SG_S242"] = "SG_2_24h_2,5mg_159"
        rename_dict["62_2d_SG_S243"] = "SG_2_48h_4,6mg_62"
        rename_dict["63_2d_SG_S244"] = "SG_2_48h_4,5mg_63"
        rename_dict["65_2d_SG_S245"] = "SG_2_48h_4,2mg_65"
        rename_dict["162_2d_SG_S246"] = "SG_2_48h_5mg_162"
        rename_dict["163_2d_SG_S264"] = "SG_2_48h_4,3mg_163"
        rename_dict["67_3d_SG_S247"] = "SG_2_72h_8,7mg_67"
        rename_dict["68_3d_SG_S265"] = "SG_2_72h_8,2mg_68"
        rename_dict["69_3d_SG_S248"] = "SG_2_72h_9,3mg_69"
        rename_dict["166_3d_SG_S249"] = "SG_2_72h_6,7mg_166"
        rename_dict["167_3d_SG_S250"] = "SG_2_72h_7,2mg_167"
        rename_dict["71_4d_SG_S266"] = "SG_2_96h_11,1mg_71"
        rename_dict["74_4d_SG_S251"] = "SG_2_96h_12,9mg_74"
        rename_dict["75_4d_SG_S252"] = "SG_2_96h_15,1mg_75"
        rename_dict["171_4d_SG_S253"] = "SG_2_96h_10,2mg_171"
        rename_dict["172_4d_SG_S254"] = "SG_2_96h_10,1mg_172"
        rename_dict["pre4_MG_S24"] = "MG_unfed_1,7mg_pre4"
        rename_dict["pre10_MG_S2"] = "MG_unfed_1,4mg_pre10"
        rename_dict["pre_8_MG_S279"] = "MG_unfed_1,6mg_pre8"
        rename_dict["4_12h_MG_S310"] = "MG_1_12h_1,8mg_4"
        rename_dict["5_12h_MG_S311"] = "MG_1_12h_1,8mg_5"
        rename_dict["104_12h_MG_S238"] = "MG_1_12h_1,8mg_104"
        rename_dict["9_24h_MG_S297"] = "MG_1_24h_1,7mg_9"
        rename_dict["10_24h_MG_S298"] = "MG_1_24h_1,6mg_10"
        rename_dict["107_24h_MG_S299"] = "MG_1_24h_1,9mg_107"
        rename_dict["13_48h_MG_S281"] = "MG_1_48h_3,2mg_13"
        rename_dict["15_48h_MG_S282"] = "MG_1_48h_3mg_15"
        rename_dict["114_48h_MG_S283"] = "MG_1_48h_2,9mg_114"
        rename_dict["20_72h_MG_S300"] = "MG_1_72h_5,5mg_20"
        rename_dict["116_72h_MG_S284"] = "MG_1_72h_6,2mg_116"
        rename_dict["118_72h_MG_S285"] = "MG_1_72h_5,9mg_118"
        rename_dict["22_96h_MG_S286"] = "MG_1_96h_10,3mg_22"
        rename_dict["24_96h_MG_S287"] = "MG_1_96h_10,2mg_24"
        rename_dict["122_96h_MG_S288"] = "MG_1_96h_9,5mg_122"
        rename_dict["155_12h_MG_S308"] = "MG_2_12h_1,8mg_155"
        rename_dict["54_12_MG_S309"] = "MG_2_12h_2,2mg_54"
        rename_dict["152_12h_MG_S276"] = "MG_2_12h_2mg_152"
        rename_dict["56_1d_MG_S277"] = "MG_2_24h_2,2mg_56"
        rename_dict["59_1d_MG_S301"] = "MG_2_24h_2,4mg_59"
        rename_dict["157_1d_MG_S278"] = "MG_2_24h_2,3mg_157"
        rename_dict["62_2d_MG_S289"] = "MG_2_48h_4,6mg_62"
        rename_dict["63_2d_MG_S290"] = "MG_2_48h_4,5mg_63"
        rename_dict["163_2d_MG_S302"] = "MG_2_48h_4,3mg_163"
        rename_dict["67_3d_MG_S291"] = "MG_2_72h_8,7mg_67"
        rename_dict["68_3d_MG_S292"] = "MG_2_72h_8,2mg_68"
        rename_dict["167_3d_MG_S293"] = "MG_2_72h_7,2mg_167"
        rename_dict["71_4d_MG_S294"] = "MG_2_96h_11,1mg_71"
        rename_dict["74_4d_MG_S295"] = "MG_2_96h_12,9mg_74"
        rename_dict["171_4d_MG_S296"] = "MG_2_96h_10,2mg_171"
        self.rename_dict = rename_dict

    def depurate(self, filter=False, nreplicates=False, matrix=False):
        if filter:
            nreplicates = int(nreplicates)
            em = Ematrix(matrix)
            with open("%s.depurated.matrix" % self.name.replace(".matrix", ""), "wt") as out:
                with open("%s.depurated_nullfiltered.matrix" % self.name.replace(".matrix", ""), "wt") as out2:
                    with open("%s.depurated_filter%iTPM.matrix" %(self.name.replace(".matrix", ""), filter), "wt") as out3:
                        out.write(self.lines[0])
                        out2.write(self.lines[0])
                        out3.write(self.lines[0])
                        for k in range(1, len(self.lines)):
                            line = self.lines[k]
                            if not line.startswith("gene") and not line.startswith("FPKM") and not line.startswith(
                                    "TPM") and not line.startswith("count"):
                                out.write(line)
                                id = line.replace("\n", "").split("\t")[0]
                                evalues = line.replace("\n", "").split("\t")[1::]
                                evalues = [float(x) for x in evalues]
                                conditions = math.floor(len(evalues)/nreplicates)

                                if not all([evalue == 0 for evalue in evalues]):
                                    out2.write(line)
                                    write = False
                                    f_evalues = [float(x) for x in em.dict[id]]
                                    for j in range(0,conditions):
                                        if all([x >= filter for x in f_evalues[j*nreplicates:(j+1)*nreplicates]]):
                                            write = True
                                    if write:
                                        out3.write(line)
        else:
            with open("%s.depurated.matrix" %self.name.replace(".matrix",""), "wt") as out:
                with open("%s.depurated_nullfiltered.matrix" % self.name.replace(".matrix", ""), "wt") as out2:
                    out.write(self.lines[0])
                    out2.write(self.lines[0])
                    for k in range(1, len(self.lines)):
                        line = self.lines[k]
                        if not line.startswith("gene") and not line.startswith("FPKM") and not line.startswith("TPM") and not line.startswith("count"):
                            out.write(line)
                            evalues = line.replace("\n","").split("\t")[1::]
                            evalues = [float(x) for x in evalues]
                            if not all([evalue == 0 for evalue in evalues]):
                                out2.write(line)

    def matrix2dict(self):
        matrix_dict = {}
        #header = ""
        for line in self.lines:
            if line.startswith("gene") == False and line.startswith("FEATURE") == False:
                ID = line.split("\t")[0]
                evalues = line.replace("\n", "").split("\t")[1::]
                matrix_dict[ID] = evalues
            #else:
                #header = header + line
        #header = header.split("\n")[0].split("\t")[1::]
        return matrix_dict

    def orderandrename(self, output, metric):
        print(self.df.head())
        order = []
        rename_df = self.df
        for key in self.rename_dict:
            order.append("%s_%s" %(self.rename_dict[key], metric))
            rename_df.columns = rename_df.columns.str.replace(key, self.rename_dict[key])

        print(rename_df.head())
        print(order)
        order_df = rename_df[order]
        order_df.to_csv(output, sep="\t")

    def rename_and_split_masigpro(self, out_label):
        header = self.lines[0]
        names = header.split("\t")
        new_names = []
        for name in names:
            print(name)
            if name != "" and "unfed" not in name:
                name = name.split("_")[0:3] + [name.split("_")[4]]
                name = "_".join(name)
                new_names.append(name)
            elif name != "" and "unfed" in name:
                name = name.split("_")[0:2] + [name.split("_")[3]]
                print(name)
                name = "_".join(name)
                new_names.append(name)
        print(new_names)
        with open("%s_full.tsv" % out_label, "wt") as out:
            out.write("\t%s\n%s" % ("\t".join(new_names), "".join(self.lines[1::])))
        with open("%s_SG.tsv" % out_label, "wt") as out:
            out.write("\t%s\n" % "\t".join(new_names)[0:55])
            for k in range(1,len(self.lines)):
                row = self.lines[k].split("\t")
                out.write("\t".join(row[0:56]))
                out.write("\n")
        with open("%s_MG.tsv" % out_label, "wt") as out:
            out.write("\t%s\n" % "\t".join(new_names)[55::])
            for k in range(1,len(self.lines)):
                row = self.lines[k].split("\t")
                out.write("%s\t%s" % (row[0], "\t".join(row[56::])))

    def split_by_tissue_masigpro(self, out_label, keep_unfed):
        header = self.lines[0].replace("\n","")
        names = [x for x in header.split("\t") if x]
        names = [x.replace('"','') for x in names]
        print(names)
        with open("%s_SG.tsv" % out_label, "wt") as out:
            if keep_unfed == True:
                out.write("\t%s\n" % "\t".join(names[0:55]))
                for k in range(1, len(self.lines)):
                    row = self.lines[k].split("\t")
                    if not all([x == "0" for x in row[1:56]]):
                        out.write("%s\t%s" % (row[0], "\t".join(row[1:56])))
                        out.write("\n")
            else:
                out.write("\t%s\n" % "\t".join(names[5:55]))
                for k in range(1,len(self.lines)):
                    row = self.lines[k].split("\t")
                    if not all([x == "0" for x in row[6:56]]):
                        out.write("%s\t%s" % (row[0], "\t".join(row[6:56])))
                        out.write("\n")
        with open("%s_MG.tsv" % out_label, "wt") as out:
            if keep_unfed == True:
                out.write("\t%s\n" % "\t".join(names[55::]))
                for k in range(1, len(self.lines)):
                    row = self.lines[k].split("\t")
                    if not all([x == "0" for x in row[56::]]):
                        out.write("%s\t%s" % (row[0], "\t".join(row[56::])))
            else:
                out.write("\t%s\n" % "\t".join(names[58::]))
                for k in range(1,len(self.lines)):
                    row = self.lines[k].split("\t")
                    if not all([x == "0" for x in row[59::]]):
                        out.write("%s\t%s" % (row[0], "\t".join(row[59::])))

    def split_by_exposure_masigpro(self, out_label, tissue, keep_unfed):
        header = self.lines[0].replace("\n", "")
        names = [x for x in header.split("\t") if x]
        names = [x.replace('"', '') for x in names]
        print(names)
        if tissue == "SG":
            with open("%s_first.tsv" % out_label, "wt") as out:
                if keep_unfed == True:
                    out.write("\t%s\n" % "\t".join(names[0:30]))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        if not all([x == "0" for x in row[1:31]]):
                            out.write("%s\t%s" % (row[0], "\t".join(row[1:31])))
                            out.write("\n")
                else:
                    out.write("\t%s\n" % "\t".join(names[5:30]))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        if not all([x == "0" for x in row[6:31]]):
                            out.write("%s\t%s" % (row[0], "\t".join(row[6:31])))
                            out.write("\n")
            with open("%s_second.tsv" % out_label, "wt") as out:
                if keep_unfed == True:
                    out.write("\t%s\t%s\n" % ("\t".join(names[0:5]), "\t".join(names[30::])))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        evalues = row[1:6]+row[31::]
                        if not all([x == "0" for x in evalues]):
                            out.write("%s\t%s" % (row[0], "\t".join(evalues)))
                else:
                    out.write("\t%s\n" % "\t".join(names[30::]))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        if not all([x == "0" for x in row[31::]]):
                            out.write("%s\t%s" % (row[0], "\t".join(row[31::])))
        elif tissue == "MG":
            with open("%s_first.tsv" % out_label, "wt") as out:
                if keep_unfed == True:
                    out.write("\t%s\n" % "\t".join(names[0:18]))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        if not all([x == "0" for x in row[1:19]]):
                            out.write("%s\t%s" % (row[0], "\t".join(row[1:19])))
                            out.write("\n")
                else:
                    out.write("\t%s\n" % "\t".join(names[3:18]))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        if not all([x == "0" for x in row[4:19]]):
                            out.write("%s\t%s" % (row[0], "\t".join(row[4:19])))
                            out.write("\n")
            with open("%s_second.tsv" % out_label, "wt") as out:
                if keep_unfed == True:
                    out.write("\t%s\t%s\n" % ("\t".join(names[0:3]), "\t".join(names[18::])))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        evalues = row[1:4]+row[19::]
                        if not all([x == "0" for x in evalues]):
                            out.write("%s\t%s" % (row[0], "\t".join(evalues)))
                else:
                    out.write("\t%s\n" % "\t".join(names[18::]))
                    for k in range(1, len(self.lines)):
                        row = self.lines[k].split("\t")
                        if not all([x == "0" for x in row[19::]]):
                            out.write("%s\t%s" % (row[0], "\t".join(row[19::])))
        else:
            print("Mistaken tissue")

    def get_experiment_design_masigpro(self, output, replicates):
        header = self.lines[0].replace("\n", "")
        names = [x for x in header.split("\t") if x]
        with open(output, "wt") as out:
            out.write("\tTime\tReplicates\tExposure.1\tExposure.2\n")
            repl = 1
            count = 0
            for name in names:
                count += 1
                if "unfed" in name:
                    time = "0"
                    exp1 = 1
                    exp2 = 1
                else:
                    time = name.split("_")[2].replace("h", "")
                    exp = name.split("_")[1]
                    if exp == "1":
                        exp1 = 1
                        exp2 = 0
                    elif exp == "2":
                        exp1 = 0
                        exp2 = 1
                out.write("%s\t%s\t%i\t%i\t%i\n" % (name, time, repl, exp1, exp2))
                if count == int(replicates):
                    count = 0
                    repl += 1

    def check_for_masigpro(self, output):
        with open(output, "wt") as out:
            out.write(self.lines[0])
            for k in range (1,len(self.lines)):
                line = self.lines[k]
                values = line.split("\t")[1::]
                values = [float(x) for x in values]
                y = np.std(values)
                if y < 0.00000000000000000000000001:
                    print(y)
                    print(k)
                    print(line.split("\t")[0])
                else:
                    out.write(self.lines[k])

    def filter(self, IDs, output):
        with open(output, "wt") as out:
            out.write(self.header)
            for key in self.dict:
                if key in IDs:
                    out.write("%s\t%s\n" %(key, "\t".join(self.dict[key])))

    def filterbyexpression(self, replicates, nconditions, threshold):
        IDs = []
        for n in range(1, len(self.lines)):
            line = self.lines[n].split("\t")
            ID_bool = True
            for j in range(0,nconditions-1):
                evalues = line[(j*replicates+1):(j*replicates+6)]
                evalues = [float(item) for item in evalues]
                mean = np.mean(evalues)
                if mean < threshold:
                    ID_bool = False
            if ID_bool == True:
                IDs.append(line[0])
        return IDs

    def filterbyexpression_MGSG(self, fasta, threshold, outlabel):
        IDs = []
        threshold = float(threshold)
        matrix_outname = "%s_FPKM.matrix" %outlabel
        fasta_outname = "%s.fasta" %outlabel
        with open(matrix_outname, "wt") as out:
            out.write(self.lines[0])
            for n in range(1, len(self.lines)):
                line = self.lines[n].split("\t")
                ID_bool = False
                for k in range(0, 11):
                    evalues = line[k * 5 + 1 : k * 5 + 6]
                    evalues = [float(item) for item in evalues]
                    if all([evalue > threshold for evalue in evalues]):
                        ID_bool = True
                for k in range(0, 11):
                    evalues = line[k * 3 + 56 : k * 3 + 59]
                    evalues = [float(item) for item in evalues]
                    if all([evalue > threshold for evalue in evalues]):
                        ID_bool = True
                if ID_bool == True:
                    out.write(self.lines[n])
                    IDs.append(line[0])

        fa = Fasta(fasta)
        fa.filter(IDs, fasta_outname)

    def filterbyexpression_sampleSheet(self, threshold, samplesheet, output):
        # Read the samplesheet
        samplesheet = pd.read_csv(samplesheet, sep="\t")

        # Map samples to groups
        sample_to_group = dict(zip(samplesheet['file'], samplesheet['group']))

        # Filter genes
        filtered_genes = []
        for gene in self.df.index:
            group_booleans = []
            group_expression_values = {}
            for sample in self.df.columns:
                if sample in sample_to_group:
                    group = sample_to_group[sample]
                    if group not in group_expression_values:
                        group_expression_values[group] = []
                    group_expression_values[group].append(self.df.loc[gene, sample])
            # Check if gene meets the condition for all groups
            print(group_expression_values)
            for values in group_expression_values.values():
                group_check = []
                for value in values:
                    group_check.append(float(value) >= float(threshold))
                group_booleans.append(all(group_check))
            if any(group_booleans):
                print("OK")
                filtered_genes.append(gene)

        # Create a new expression matrix with filtered genes
        filtered_expression_df = self.df.loc[filtered_genes, :]

        # Save the filtered expression matrix to a new file
        filtered_expression_df.to_csv(output, sep='\t')

    def genestoqPCR(self, cexcel, output):
        genes_5FPKM = {}
        cds_dict = {}
        transcript_dict = {}
        cexcel = CompleteExcel(cexcel)
        for n in range(1,len(self.lines)):
            ID = self.lines[n].split("\t")[0]
            evalues = self.lines[n].replace("\n","").split("\t")[1::]
            if all(float(x) > 5 for x in evalues):
                genes_5FPKM[ID] = evalues
        df = pd.DataFrame.from_dict(genes_5FPKM, orient="index", columns=self.df.columns)
        df = df.astype(float)
        df['mean'] = df.mean(axis=1)
        df = df.sort_values("mean", ascending=False)
        df_head = df.head(n=30)
        df_head.to_csv("%s.matrix" %output, sep="\t")
        IDs = df_head.index
        for n in range(1,len(cexcel.lines)):
            line = CompleteExcel.Line(cexcel.lines[n])
            if line.transcriptID in IDs:
                transcript_dict[line.transcriptID] = line.transcript_seq
                cds_dict[line.transcriptID] = line.cds_seq
        with open("%s.fa" % output, "wt") as out_fa:
            with open("%s.cds" % output, "wt") as out_cds:
                for ID in IDs:
                    out_fa.write(">%s\n%s\n" % (ID, transcript_dict[ID]))
                    out_cds.write(">%s\n%s\n" % (ID, cds_dict[ID]))

    def hierarchical_clustering(self, out_label, method):
        df_hier = self.df.transpose()
        df_hier.index.name = None
        methods = ["ward", "single", "complete", "average", "weighted", "centroid", "median"]
        if method in methods:
            plt.figure(figsize=(10, 5))
            plt.title("Dendrogram %s" %method)
            link = shc.linkage(df_hier, method=method)
            dend = shc.dendrogram(link, labels=df_hier.index, leaf_rotation=90)
            plt.tight_layout()
            plt.savefig("%s_%s_dendrogram.tiff" %out_label %(out_label, method), dpi=300)
            plt.close()
        elif method == "all":
            methods = ["ward", "single", "complete", "average", "weighted", "centroid", "median"]
            for meth in methods:
                plt.figure(figsize=(10, 5))
                plt.title("Dendrogram %s" % meth)
                link = shc.linkage(df_hier, method=meth)
                dend = shc.dendrogram(link, labels=df_hier.index, leaf_rotation=90)
                plt.tight_layout()
                plt.savefig("%s_%s_dendrogram.tiff" %(out_label, meth), dpi=300)
                plt.close()
        else:
            print("Please select a valid method for hierarchical clustering")

    def correlation_matrix(self, out_label):
        corr = self.df.corr()
        corr.to_csv("%s_corr.csv" % out_label)
        plt.figure(figsize=(7, 7))

        #ax = sns.set(font_scale=0.5)

        ax = sns.heatmap(
            corr,
            vmin=-1, vmax=1, center=0,
            cmap=sns.diverging_palette(20, 220, n=200),
            square=True
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            horizontalalignment='right'
        )
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig("%s_corr.tiff" %out_label, dpi=300)
        plt.close()

    def get_vectors(self, sep, cn, fn):
        color_vector = []
        form_vector = []
        for index, row in self.df.transpose().iterrows():
            ID = index.split(sep)
            color_vector.append(ID[cn])
            form_vector.append(ID[fn])
        return color_vector, form_vector

    def PCA(self, out_label, color_vector, form_vector):
        df_transposed_standarized = StandardScaler().fit_transform(self.df.transpose())
        PCA_data = PCA(n_components=2).fit_transform(df_transposed_standarized)
        df_PCA = pd.DataFrame(data=PCA_data, index=self.df.transpose().index, columns=["component1", "component2"])
        df_PCA["group"] = color_vector
        df_PCA["tissue"] = form_vector
        df_PCA = df_PCA.sort_values(by=["group"], ascending=False)

        fig, ax = plt.subplots(facecolor='w')

        plt.figure(figsize=(7, 7))

        ax = sns.set(font_scale=1)
        ax = sns.scatterplot(data=df_PCA, x="component1", y="component2", hue="group", style="tissue",
                             palette="coolwarm",
                             edgecolor="black")
        for key, row in df_PCA.iterrows():
            ax.annotate(key, xy=(row['component1'], row['component2']))
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig("%s_PCA_withlabels.tiff" %out_label, dpi=300)
        plt.close()

        plt.figure(figsize=(7, 7))
        ax = sns.set(font_scale=1)
        ax = sns.scatterplot(data=df_PCA, x="component1", y="component2", hue="group", style="tissue",
                             palette="coolwarm",
                             edgecolor="black")
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig("%s_PCA.tiff" %out_label, dpi=300)
        plt.close()

    def MDS(self, out_label, color_vector, form_vector):
        df_transposed_standarized = StandardScaler().fit_transform(self.df.transpose())
        mds_data = MDS(metric=True, random_state=0).fit_transform(df_transposed_standarized)
        df_mds = pd.DataFrame(data=mds_data, index=self.df.transpose().index, columns=["X", "Y"])
        df_mds["group"] = color_vector
        df_mds["tissue"] = form_vector
        df_mds = df_mds.sort_values(by=["group"], ascending=False)

        fig, ax = plt.subplots(facecolor='w')

        plt.figure(figsize=(7, 7))
        ax = sns.set(font_scale=1)
        ax = sns.scatterplot(data=df_mds, x="X", y="Y", hue="group", style="tissue",
                             palette="coolwarm",
                             edgecolor="black")
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig("%s_MDS.tiff" %out_label, dpi=300)
        plt.close()

        nmds_data = MDS(metric=False, random_state=0).fit_transform(df_transposed_standarized)
        df_nmds = pd.DataFrame(data=nmds_data, index=self.df.transpose().index, columns=["X", "Y"])
        df_nmds["group"] = color_vector
        df_nmds["tissue"] = form_vector
        df_nmds = df_nmds.sort_values(by=["group"], ascending=False)

        fig, ax = plt.subplots(facecolor='w')

        plt.figure(figsize=(7, 7))
        ax = sns.set(font_scale=1)
        ax = sns.scatterplot(data=df_nmds, x="X", y="Y", hue="group", style="tissue",
                             palette="coolwarm",
                             edgecolor="black")
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig("%s_nMDS.tiff" %out_label, dpi=300)
        plt.close()

    def expression_plot_by_gene(self, output_folder, label):
        header = self.header.split("\t")
        for gene in self.dict:
            D = dict()
            for k in range(0,len(self.dict[gene])):
                D[header[k+1]] = float(self.dict[gene][k])
            out = "%s/%s.tiff" %(output_folder, gene)
            barplot_from_dict(D,out,gene, label)

    def expression_plot_by_gene_error(self, error_tsv, output_folder, label):
        header = self.header.split("\t")
        errormat = Ematrix(error_tsv)
        for gene in self.dict:
            D = dict()
            D_error = dict()
            for k in range(0,len(self.dict[gene])):
                D[header[k+1]] = float(self.dict[gene][k])
                D_error[header[k + 1]] = float(errormat.dict[gene][k])
            out = "%s/%s.png" %(output_folder, gene)
            barplot_from_dict_error(D, D_error, out,gene, label)

    def remove_columns_containing_string(self, string, output):
        df = self.df[self.df.columns.drop(list(self.df.filter(regex=string)))]
        df.to_csv(output, sep="\t")

    def get_ranking(self, out_label, n=100):
        evalues = []
        genes = []
        for k in range (1, len(self.lines)):
            gene = self.lines[k].split("\t")[0]
            genes.append(gene)
            line_evalues = self.lines[k].replace("\n","").split("\t")[1::]
            line_evalues = [float(x) for x in line_evalues]
            evalues.append(line_evalues)

        evalues = np.array(evalues)

        print(evalues)
        rankings =scipy.stats.rankdata(-evalues, axis=-0, method="ordinal")
        print(rankings)
        with open("%s_rankings.tsv" %out_label, "wt") as out:
            out.write(self.header)
            for k in range(0, len(rankings)):
                out.write("%s\t%s\n" %(genes[k], "\t".join([str(x) for x in rankings[k]])))

        rankings_matrix = np.array(rankings, dtype=float)

        # Get the top 100 genes for each sample
        top_genes_per_sample = []

        for i in range(0, len(self.header.split("\t")) - 1):  # Iterate over samples
            sample_rankings = rankings_matrix[:, i]
            top_genes_indices = np.argwhere(sample_rankings < int(n)).flatten()  # Get indices of the top 100 genes
            top_genes_per_sample.append(top_genes_indices)

        print(top_genes_per_sample)
        genes_dict = dict()
        k = 0
        for sample in top_genes_per_sample:
            k = k + 1
            sam = self.header.replace("\n","").split("\t")[k]
            genes_dict[sam] = []
            for gene_index in sample:
                genes_dict[sam].append(genes[int(gene_index)])

        df = pd.DataFrame.from_dict(genes_dict)

        # Save the DataFrame to a CSV file
        df.to_csv("%s_top%igenes.csv" %(out_label, int(n)), index=False)

class DE:
    def __init__(self, file, excel_transcriptIDs=False):
        with open(file) as de:
            self.lines = list(de)
        self.IDsUPfc = {}
        self.IDsDOWNfc = {}
        self.IDsUPpvalue = {}
        self.IDsDOWNpvalue = {}
        for n in range(1, len(self.lines)):
            line = DE.Line(self.lines[n])
            if excel_transcriptIDs == False:
                if line.fc >= 0:
                    self.IDsUPfc[line.ID] = "%.4f" %line.fc
                    self.IDsUPpvalue[line.ID] = "%.4f" %line.pvalue
                else:
                    self.IDsDOWNfc[line.ID] = "%.4f" %line.fc
                    self.IDsDOWNpvalue[line.ID] = "%.4f" %line.pvalue
            else:
                if line.ID in excel_transcriptIDs:
                    if line.fc >= 0:
                        self.IDsUPfc[line.ID] = "%.4f" %line.fc
                        self.IDsUPpvalue[line.ID] = "%.4f" %line.pvalue
                    else:
                        self.IDsDOWNfc[line.ID] = "%.4f" %line.fc
                        self.IDsDOWNpvalue[line.ID] = "%.4f" %line.pvalue

        self.IDsUP = list(self.IDsUPfc.keys())
        self.IDsDOWN = list(self.IDsDOWNfc.keys())
        self.IDs = self.IDsUP.extend(self.IDsDOWN)

    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.ID = self.row[0].replace('"','')
            self.pvalue = float(self.row[2])
            self.fc = np.log2(float(self.row[3]))

class Excel:
    def __init__(self, file):
        with open(file) as ex:
            self.lines = list(ex)
        self.IDs = self.getIDlist()
        self.dict = self.excel2dict()

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.cdsID = self.row[0]
            self.seq = self.row[1]
            self.hitID = gethitID(self.row)
            self.identity = self.row[3]
            self.coverage = self.row[4]
            self.bitscore = self.row[5]
            self.description = self.row[6]
            self.db = self.row[7]
            self.ex_count = self.row[8:96]
            self.fpkm = self.row[96:184]

    def getnames(self):
        names = []
        for n in range(1, len(self.lines)):
            line = Excel.Line(self.lines[n])
            names.append(line.hitID)
        names = deduplicate_list(names)
        return names

    def getIDlist(self):
        IDlist = []
        for n in range(1, len(self.lines)):
            line = Excel.Line(self.lines[n])
            IDlist.append(line.cdsID)
        return IDlist

    def excel2dict(self):
        ex_dict = {}
        for n in range(1, len(self.lines)):
            line = Excel.Line(self.lines[n])
            ex_dict[line.cdsID] = self.lines[n]
        return ex_dict

class Uniprot:
    def __init__(self, file, obo_dict):
        with open(file) as uni:
            self.lines = list(uni)
        self.goterm_dict, self.goterm_des_dict, self.keyword_dict, self.keyword_des_dict, self.family_dict = self.uniprot2dict(obo_dict)

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.hitID = self.row[-1]
            self.entry = self.row[0]
            self.entry_name = self.row[1]
            self.protein_name = self.row[2]
            self.organism = self.row[3]
            self.GObp = self.row[4]
            self.GOcc = self.row[5]
            self.GO = self.row[6]
            self.GOmf = self.row[7]
            self.GOterm = self.row[8]
            self.keyword = self.row[9]
            self.keywordID = self.row[10]
            self.family = self.row[11]


    def uniprot2dict(self,obo_dict):
        goterm_dict = {}
        goterm_des_dict = {}
        kw_dict = {}
        kw_des_dict = {}
        family_dict = {}
        for n in range(1,len(self.lines)):
            uniprot_line = Uniprot.Line(self.lines[n])
            if uniprot_line.GOterm != "":
                goterm_dict[uniprot_line.hitID] = uniprot_line.GOterm
                goterm_des_dict[uniprot_line.hitID] = "; ".join([obo_dict[x] for x in uniprot_line.GOterm.split("; ")])
            else:
                goterm_dict[uniprot_line.hitID] = "None"
                goterm_des_dict[uniprot_line.hitID] = "None"
            if uniprot_line.keywordID != "":
                kw_dict[uniprot_line.hitID] = uniprot_line.keywordID
                kw_des_dict[uniprot_line.hitID] = uniprot_line.keyword
            else:
                kw_dict[uniprot_line.hitID] = "None"
                kw_des_dict[uniprot_line.hitID] = "None"
            if uniprot_line.family != "":
                family_dict[uniprot_line.hitID] = uniprot_line.family

        return goterm_dict, goterm_des_dict, kw_dict, kw_des_dict, family_dict

    def getAnnotated(self):
        hitID_list = []
        for n in range(1, len(self.lines)):
            line = Uniprot.Line(self.lines[n])
            if line.GOterm != "":
                hitID_list.append(line.hitID)
        return hitID_list

class Blast2GO:
    def __init__(self, file, obo_dict):
        with open(file) as bg:
            self.lines = list(bg)
        self.dict, self.des_dict = self.blast2GO2dict(obo_dict)

    def blast2GO2dict(self, obo_dict):
        goterm_dict = {}
        goterm_des_dict = {}
        for n in range(1, len(self.lines)):
            ID = self.lines[n].split("\t")[2]
            goterms = self.lines[n].split("\t")[9].replace("C:", "").replace("P:", "").replace("F:", "")
            if goterms != "" and goterms != "no GO terms":
                goterm_dict[ID] = goterms
                goterm_des_dict[ID] = "; ".join([obo_dict[x] for x in goterms.split("; ")])
            else:
                goterm_dict[ID] = "None"
                goterm_des_dict[ID] = "None"
        return goterm_dict, goterm_des_dict

class Fasta:
    def __init__(self, file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = self.fasta2dict()
        self.full_dict = self.fasta2fulldict()
        self.IDs = self.dict.keys()
        self.len_dict = self.get_len_dict()

    def split(self, n, output):
        file_number = 1
        count = 0
        out = open("%s_%s.fa" % (output, file_number), "wt")
        for key in self.dict:
            out.write(">%s\n%s\n" % (key, self.dict[key]))
            count += 1
            if count >= int(n):
                out.close()
                file_number += 1
                count = 0
                out = open("%s_%s.fa" % (output, file_number), "wt")

    def fasta2dict(self):
        fasta_dict = dict()
        seqid = ""
        for i in range(0, len(self.lines)):
            if self.lines[i].startswith('>'):
                if seqid != "":
                    fasta_dict[seqid] = seq
                seqid = self.lines[i].split(' ')[0].replace(">", "").replace("\n", "")
                seq = ""
            else:
                seq += self.lines[i].replace("\n", "")
        if seqid != "":
            fasta_dict[seqid] = seq
        return fasta_dict

    def fasta2fulldict(self):
        fasta_dict = dict()
        seqid = ""
        for i in range(0, len(self.lines)):
            if self.lines[i].startswith('>'):
                if seqid != "":
                    fasta_dict[seqid] = seq
                seqid = self.lines[i].replace("\n", "").replace(">","")
                seq = ""
            else:
                seq += self.lines[i].replace("\n", "")
        if seqid != "":
            fasta_dict[seqid] = seq
        return fasta_dict

    def filterlongestisoform(self,output):
        form_dict = {}
        for key in self.dict:
            form_ID = "_".join(key.split("_")[0:-1])
            if form_ID not in form_dict:
                form_dict[form_ID] = key
            else:
                if len(self.dict[key]) > len(self.dict[form_dict[form_ID]]):
                    form_dict[form_ID] = key
        with open(output,"wt") as out:
            for key in form_dict:
                out.write(">%s\n%s\n" %(form_dict[key], self.dict[form_dict[key]]))

    def filter(self,IDlist,output):
        IDlist_2 = []
        for ID in IDlist:
            # if "." in ID:
            #      IDlist_2.append(ID.split(".")[0])
            # else:
            #      IDlist_2.append(ID)
            IDlist_2.append(ID)

        with open(output, "wt") as out:
            for key in self.dict:
                # if "." in key:
                #      id = key.split(".")[0]
                # else:
                #      id = key
                id = key
                if id in IDlist or key in IDlist_2 or id in IDlist_2:
                    out.write(">%s\n%s\n" %(id, self.dict[key]))

    def len_distribution(self,out_label):
        len_list = []

        for key in self.dict:
            len_list.append(len(self.dict[key]))


        max = np.max(len_list)
        p25 = np.percentile(len_list, 25)
        p95 = np.percentile(len_list, 95)
        p75 = np.percentile(len_list, 75)
        p50 = np.percentile(len_list, 50)
        p5 = np.percentile(len_list, 5)
        min = np.min(len_list)
        mean = np.mean(len_list)
        median = np.median(len_list)
        n50 = calculate_N50(len_list)

        with open("%s_lengthstats.txt" %out_label, "wt") as out:
            out.write("max\tmin\tmean\tmedian\tn50\tp95\tp75\tp50\tp25\tp5\n")
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (max, min, mean, median, n50, p95, p75, p50, p25, p5))

        fig = plt.figure()
        sns.set_theme()
        plt.hist(len_list, bins=50)
        #plt.axvline(x=n50, color='green', linestyle='--', label="N50")
        plt.axvline(x=mean, color='red', linestyle='--', label="Mean")
        plt.legend(loc="upper right")
        fig.savefig("%s_displot.png" %out_label)

        fig = plt.figure()
        sns.set_theme()
        plt.hist(len_list, bins=50, log=True)
        #plt.axvline(x=n50, color='green', linestyle='--', label="N50")
        plt.axvline(x=mean, color='red', linestyle='--', label="Mean")
        plt.legend(loc="upper right")
        fig.savefig("%s_displot_log.png" % out_label)

    def get_len_dict(self):
        len_dict = dict()
        for key in self.dict:
            len_dict[key] = len(self.dict[key])

        return len_dict

    def filter_by_length(self, output, max_length=1000000000000000, min_length=0):
        with open(output, "wt") as out:
            for key in self.dict:
                if len(self.dict[key]) >= min_length and len(self.dict[key]) <= max_length:
                    out.write(">%s\n%s\n" %(key, self.dict[key]))

    def filter_by_name(self,name,output):
        with open(output, "wt") as out:
            for key in self.full_dict:
                if name in key.lower():
                    out.write(">%s\n%s\n" %(key, self.full_dict[key]))

    def deduplicate(self, output):
        with open(output, "wt") as out:
            for key in self.full_dict:
                out.write(">%s\n%s\n" %(key, self.full_dict[key]))

    def sort_apollo_fasta(self, output):
        names = []
        chroms = []
        starts = []
        ends = []
        coords = []
        domains = []
        for key in self.full_dict:
            name = key.split(" ")[0]
            print(name)
            chrom = int(key.split("[Iricinus_Scaffold_")[1].split(":")[0])
            start = int(key.split("[")[1].split("strand]")[0].split(":")[1].split("-")[0])
            end = int(key.split("[")[1].split("strand]")[0].split(":")[1].split(" ")[0].split("-")[1])
            coord = "Iricinus_Scaffold_%i:%i-%i" %(chrom, start, end)
            domain = int(key.split("K=")[1])
            names.append(name)
            chroms.append(chrom)
            starts.append(start)
            ends.append(end)
            coords.append(coord)
            domains.append(domain)
        df = pd.DataFrame.from_dict({"ID":names, "Chrom":chroms, "Start":starts, "End":ends, "Coords":coords, "K domains":domains})
        df_sorted = df.sort_values(['Chrom', 'Start', "End"])
        df_sorted.to_csv(output, sep="\t")

    def fasta2tsv(self, output):
        with open(output, "wt") as out:
            for key in self.dict:
                out.write("%s\t%s\n" %(key, self.dict[key]))

    def trim(self, sep, output):
        with open(output,"wt") as out:
            for key in self.full_dict:
                if sep in key:
                    name= key.split(sep)[0]
                else:
                    name=key
                out.write(">%s\n%s\n" %(name, self.full_dict[key]))

    def intersect(self, fa2, output, reverse=False):
        fa = Fasta(fa2)

        ids = []
        for id in fa.dict.keys():
            if "." in id:
                id = id.split(".")[0]
            ids.append(id)
        with open(output,"wt") as out:
            for key in self.full_dict:
                id = key.split(" ")[0]
                if "." in id:
                    id = id.split(".")[0]
                if reverse == False:
                    if id in ids:
                        out.write(">%s\n%s\n" %(key, self.full_dict[key]))
                else:
                    if id not in ids:
                        out.write(">%s\n%s\n" %(key, self.full_dict[key]))

    def print_lengths(self):
        for key in self.len_dict:
            print("%s: %i" %(key, self.len_dict[key]))

class Interpro:
    def __init__(self, file, obo_dict):
        with open(file) as ip:
            self.lines = list(ip)
        self.ipro_dict, self.ipro_des_dict, self.goterm_dict, self.goterm_des_dict = self.interpro2dict(obo_dict)

    class Line:
        def __init__(self, line):
            fields = line.replace("\n", "").split("\t")
            self.ID = fields[0]
            self.Analysis = fields[3]
            self.sign_acc = fields[4]
            self.sign_des = fields[5]
            self.ipro_annotation = fields[11]
            self.ipro_description = fields[12]
            if len(fields) > 13:
                if fields[13] != "-":
                    self.goterm = fields[13].replace("|","; ")
                else:
                    self.goterm = "None"
            else:
                self.goterm = "None"

    def interpro2dict(self, obo_dict):
        ipro_dict = {}
        ipro_des_dict = {}
        goterm_dict = {}
        goterm_des_dict = {}
        for line in self.lines:
            ipro_line = Interpro.Line(line)
            if ipro_line.ID not in ipro_dict and ipro_line.ipro_annotation != "-":
                ipro_dict[ipro_line.ID] = [ipro_line.ipro_annotation]
                ipro_des_dict[ipro_line.ID] = [ipro_line.ipro_description]
            elif ipro_line.ID in ipro_dict and ipro_line.ipro_annotation != "-":
                if ipro_line.ipro_annotation not in ipro_dict[ipro_line.ID]:
                    ipro_dict[ipro_line.ID].append(ipro_line.ipro_annotation)
                    ipro_des_dict[ipro_line.ID].append(ipro_line.ipro_description)
            if ipro_line.ID not in goterm_dict and ipro_line.goterm != "None":
                goterm_dict[ipro_line.ID] = [ipro_line.goterm]
                goterm_des_dict[ipro_line.ID] = [obo_dict[x] for x in ipro_line.goterm.split("; ")]
            elif ipro_line.ID in goterm_dict and ipro_line.goterm != "None":
                if ipro_line.goterm not in goterm_dict[ipro_line.ID]:
                    goterm_dict[ipro_line.ID].append(ipro_line.goterm)
                    goterm_des_dict[ipro_line.ID] += [obo_dict[x] for x in ipro_line.goterm.split("; ")]

        for key in ipro_dict:
            ipro_dict[key] = "; ".join(deduplicate_list("; ".join(ipro_dict[key]).split("; ")))
            ipro_des_dict[key] = "; ".join(deduplicate_list("; ".join(ipro_des_dict[key]).split("; ")))
        for key in goterm_dict:
            goterm_dict[key] = "; ".join(deduplicate_list("; ".join(goterm_dict[key]).split("; ")))
            goterm_des_dict[key] = "; ".join(deduplicate_list("; ".join(goterm_des_dict[key]).split("; ")))

        return ipro_dict, ipro_des_dict, goterm_dict, goterm_des_dict

class SignalP:
    def __init__(self, file):
        with open(file) as sp:
            self.lines = list(sp)
        self.header = self.lines[0:2]
        self.dict = self.signalP2dict()

    def signalP2dict(self):
        sp_dict = {}
        for n in range(2, len(self.lines)):
            ID = self.lines[n].split("\t")[0]
            signal = self.lines[n].split("\t")[1]
            if signal == "OTHER":
                signal = "None"
            sp_dict[ID] = signal
        return sp_dict

class SignalP4:
    def __init__(self, file, tmhmm_mature):
        with open(file) as sp:
            self.lines = list(sp)
        self.dict = self.signalP42dict()

    def signalP42dict(self):
        sp_dict = {}
        count = 0
        for n in range(2,len(self.lines)):
            ID = self.lines[n].split("  ")[0]
            sp = self.lines[n].split("  ")[8][-1]
            if sp == "Y":
                count += 1
                sp = "Yes"
            else:
                sp = "None"
            sp_dict[ID] = sp
        print("%s total of signal peptides" %count)
        return sp_dict

class Phobius:
    def __init__(self, file):
        with open(file) as ph:
            self.lines = list(ph)
        self.dict = self.getSecreted()

    class Line:
        def __init__(self,line):
            self.row = "\t".join(line.split()).split("\t")
            self.ID = self.row[0]
            self.tmhmm = self.row[1]
            self.sp = self.row[2]
            self.prediction = self.row[3]

    def getSecreted(self):
        secreted_dict = {}
        count = 0
        for n in range(1,len(self.lines)):
            line = Phobius.Line(self.lines[n])
            if line.tmhmm == "0" and line.sp == "Y":
                count += 1
                secreted_dict[line.ID] = "Yes"
            else:
                secreted_dict[line.ID] = "None"
        perc = count/(len(self.lines) - 1) * 100
        #print("%i Phobius secreted proteins; %f of the total of %i" %(count,perc,len(self.lines)-1))
        return secreted_dict

class TMHMM:
    def __init__(self, file):
        with open(file) as tm:
            self.lines = list(tm)
        self.TMHlines = self.getTMHlines()
        self.dict = self.tmh2dict()

    def getTMHlines(self):
        TMHlines = []
        for line in self.lines:
            if "Number of predicted TMHs" in line:
                TMHlines.append(line)
        return TMHlines

    def tmh2dict(self):
        tmh_dict = {}
        for line in self.TMHlines:
            ID = line.split(" Number of predicted TMHs:  ")[0].replace("# ", "")
            tmh_n = line.split(" Number of predicted TMHs:  ")[-1].replace("\n", "")
            if tmh_n == "0":
                tmh_n = "None"
            tmh_dict[ID] = tmh_n
        return tmh_dict

class GOobo:
    def __init__(self, file):
        with open(file) as obo:
            self.lines = list(obo)
        self.dict = {}
        self.namespace_dict = {}
        for n in range(0,len(self.lines)):
            if "[Term]" in self.lines[n]:
                ID = self.lines[n+1].split(": ")[1].replace("\n","")
                name = self.lines[n+2].split(": ")[1].replace("\n","")
                namespace = self.lines[n+3].split(": ")[1].replace("\n","")
                self.dict[ID] = name
                self.namespace_dict[ID] = namespace
            elif "alt_id:" in self.lines[n]:
                ID = self.lines[n].split(": ")[1].replace("\n","")
                self.dict[ID] = name
                self.namespace_dict[ID] = namespace

class KEGobo:
    def __init__(self,file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = dict()
        for line in self.lines:
            ID = line.split("\t")[0].split(":")[1]
            description = line.split("\t")[1].replace("\n","")
            self.dict[ID] = description

class EnrichmentOut:
    def __init__(self, file):
        with open(file) as enr:
            self.lines = list(enr)

    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.GO = self.row[0]
            self.description = self.row[1]
            self.inbackground = int(self.row[2])
            self.notinbackground = int(self.row[3])
            self.percentageinbackground = self.inbackground / (self.inbackground + self.notinbackground) * 100
            self.ingenelist = int(self.row[4])
            self.notingenelist = int(self.row[5])
            self.percentageingenelist = self.ingenelist / (self.ingenelist + self.notingenelist) * 100
            self.pvalue = float(self.row[6])
            self.fdr = float(self.row[7])
            self.oddsratio = float(self.row[8])
            self.RE = float(self.row[9])
            self.type = self.row[10]

    def parsedescription(self, des_dict, output):
        with open(output,"wt") as out:
            out.write(self.lines[0])
            for n in range(1,len(self.lines)):
                line = EnrichmentOut.Line(self.lines[n])
                des = des_dict[line.ID]
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line.ID,des,line.inbackground,line.ingenelist,line.notingenelist,line.pvalue,line.fdr,line.re,line.type,line.genes))

    def getGO(self):
        GO_dict = {}
        for n in range(1,len(self.lines)):
                line = EnrichmentOut.Line(self.lines[n])
                GO_dict[line.ID] = [line.des, line.type]
        return GO_dict
    
    def depurate(self, GO_description_file, output):
        go_descriptions = GOobo(GO_description_file)

        with open(output, "wt") as out:
            out.write(self.lines[0])
            for j in range(1, len(self.lines)):
                en_line = EnrichmentOut.Line(self.lines[j])
                if int(en_line.ingenelist) > 0:
                    en_line.row[1] = go_descriptions.dict[en_line.GO]
                    out.write("\t".join(en_line.row))

    def depurateKEG(self, KEG_description_file, output):
        go_descriptions = KEGobo(KEG_description_file)

        with open(output, "wt") as out:
            out.write(self.lines[0])
            for j in range(1, len(self.lines)):
                en_line = EnrichmentOut.Line(self.lines[j])
                if int(en_line.ingenelist) > 0:
                    en_line.row[1] = go_descriptions.dict[en_line.GO]
                    out.write("\t".join(en_line.row))

    def getReviGOFile(self, output):
        with open(output, "wt") as out:
            for j in range(1, len(self.lines)):
                en_line = EnrichmentOut.Line(self.lines[j])
                out.write("%s\t%s\n" %(en_line.ID, en_line.pvalue))

    def barplot_rawdata(self,outname, number):
        bar_dict = dict()
        for j in range(1, len(self.lines)):
            en_line = EnrichmentOut.Line(self.lines[j])
            bar_dict[en_line.description.split(",")[0]] = int(en_line.ingenelist)
        bar_dict = order_dictionary(bar_dict)
        barhplot_from_dict(bar_dict,outname,"Absolute value", number)

    def barplot_proportion(self,outname, number, filter_n=False):
        bar_dict = dict()
        for j in range(1, len(self.lines)):
            en_line = EnrichmentOut.Line(self.lines[j])
            if filter_n and float(en_line.ingenelist) >= filter_n:
                bar_dict[en_line.description.split(",")[0]]=float(en_line.ingenelist)/float(en_line.inbackground) * 100
            elif filter_n==False:
                bar_dict[en_line.description.split(",")[0]] = float(en_line.ingenelist) / float(
                    en_line.inbackground) * 100
        bar_dict = order_dictionary(bar_dict)
        barhplot_from_dict(bar_dict, outname,"inGeneList/inBackground (%)", number)

    def barhplot_bujun(D, out, number):
        transcript_list = []
        ingene_list = []
        background_list = []

        if len(D) > number:
            n = 0
            for key in D:
                transcript_list.append(key)
                ingene_list.append(D[key][0])
                background_list.append(D[key][1])
                n += 1
                if n == int(number):
                    break
        else:
            for key in D:
                transcript_list.append(key)
                ingene_list.append(D[key][0])
                background_list.append(D[key][1])

        df = pd.DataFrame({'% in gene list': ingene_list, '% in background': background_list}, index=transcript_list)
        sns.set_theme()
        ax = df.plot.barh()
        ax.set_ylabel("GO term")
        fig = ax.get_figure()
        fig.set_figheight(10)
        fig.set_figwidth(20)
        plt.tight_layout()
        fig.savefig(out)
        plt.close()

    def parse2GOfigure(self, output):
        with open(output,"wt") as out:
            out.write("GO term\tenrichment_P-value\tOdds ratio\n")
            for k in range(1, len(self.lines)):
                line = EnrichmentOut.Line(self.lines[k])
                if line.type == "enriched" and line.pvalue < 0.05:
                    out.write("%s\t%s\t%s\n" %(line.GO, line.pvalue, line.oddsratio))

class ConditionCV:
    def __init__(self,file):
        with open(file) as f:
            self.lines = list(f)
        self.conditions = self.lines[0].split("\t")[1:23]
        self.df = pd.read_csv(file, sep="\t")

    def expressionfilter(self,n):
        conditions_dict = {}
        print(self.lines[0].split("\t")[23:28])
        print(self.lines[0].split("\t")[11*3+45:11*3+48])
        for j in range(0,22):
            conditions_dict[self.conditions[j]] = []
        for k in range(1,len(self.lines)):
            ID = self.lines[k].split("\t")[0].split(".")[0]
            for j in range(0,11):
                values = self.lines[k].split("\t")[j*5+23:j*5+28]
                if all([float(value) >= n for value in values]):
                    conditions_dict[self.conditions[j]].append(ID)
            for j in range(11,22):
                values = self.lines[k].split("\t")[j*3+45:j*3+48]
                if all([float(value) >= n for value in values]):
                    conditions_dict[self.conditions[j]].append(ID)
        return conditions_dict




    def getEnrichmentInput(self, output_label, n, reference_label=""):
        list_of_files=[]
        for condition in self.conditions:
            IDs = []
            df = self.df.sort_values(condition, ascending=False)
            label = "%s_%s" %(output_label,"_".join(condition.split()))
            name = "%s_mostCV.txt" %(label)
            list_of_files.append(name)
            with open(name, "wt") as out:
                for j in range(0, n):
                    ID = df.iloc[j]["ID"].split(".")[0]
                    IDs.append(ID)
                    out.write("%s\n" %ID)
                df = self.df.sort_values(condition)
                name = "%s_leastCV.txt" %(label)
                list_of_files.append(name)
            with open(name, "wt") as out:
                for j in range(0, n):
                    ID = df.iloc[j]["ID"].split(".")[0]
                    out.write("%s\n" %ID)
                    IDs.append(ID)
            if reference_label != "":
                intersect_background(IDs,"%s_GO_all.txt" %reference_label,"%s_GO_all.txt" %label)
                intersect_background(IDs,"%s_KW_all.txt" %reference_label,"%s_KW_all.txt" %label)
                intersect_background(IDs,"%s_IP_all.txt" %reference_label,"%s_IP_all.txt" %label)
                intersect_background(IDs,"%s_family_all.txt" %reference_label,"%s_family_all.txt" %label)
        return list_of_files

    def getEnrichmentInputAnnotated(self, output_label, n, reference_label, completeExcel, pwd, expressionfiltern=False):
        cExcel = CompleteExcel(completeExcel)
        os.system("mkdir %s/GO %s/KW %s/IP %s/family" % (pwd, pwd, pwd, pwd))
        if expressionfiltern != False:
            conditions_dict = self.expressionfilter(expressionfiltern)
        for condition in self.conditions:
            IDs_GO = []
            IDs_KW = []
            IDs_IP = []
            IDs_family = []
            label = "%s_%s" %(output_label,"_".join(condition.split()))
            df = self.df.sort_values(condition, ascending=False)
            for j in range(0,df.shape[0]):
                ID = df.iloc[j]["ID"].split(".")[0]
                if expressionfiltern == False:
                    if ID in cExcel.IDs_GO and len(IDs_GO) < n:
                        IDs_GO.append(ID)
                    if ID in cExcel.IDs_KW and len(IDs_KW) < n:
                        IDs_KW.append(ID)
                    if ID in cExcel.IDs_IP and len(IDs_IP) < n:
                        IDs_IP.append(ID)
                    if ID in cExcel.IDs_family and len(IDs_family) < n:
                        IDs_family.append(ID)
                else:
                    if ID in cExcel.IDs_GO and len(IDs_GO) < n and ID in conditions_dict[condition]:
                        IDs_GO.append(ID)
                    if ID in cExcel.IDs_KW and len(IDs_KW) < n and ID in conditions_dict[condition]:
                        IDs_KW.append(ID)
                    if ID in cExcel.IDs_IP and len(IDs_IP) < n and ID in conditions_dict[condition]:
                        IDs_IP.append(ID)
                    if ID in cExcel.IDs_family and len(IDs_family) < n and ID in conditions_dict[condition]:
                        IDs_family.append(ID)
            with open("%s_GO_mostCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_GO))
            with open("%s_KW_mostCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_KW))
            with open("%s_IP_mostCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_IP))
            with open("%s_family_mostCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_family))
            IDs_GO_global = IDs_GO
            IDs_KW_global = IDs_KW
            IDs_IP_global = IDs_IP
            IDs_family_global = IDs_family
            IDs_GO = []
            IDs_KW = []
            IDs_IP = []
            IDs_family = []
            df = self.df.sort_values(condition)

            for j in range(0,df.shape[0]):
                ID = df.iloc[j]["ID"].split(".")[0]
                if expressionfiltern == False:
                    if ID in cExcel.IDs_GO and len(IDs_GO) < n:
                        IDs_GO.append(ID)
                    if ID in cExcel.IDs_KW and len(IDs_KW) < n:
                        IDs_KW.append(ID)
                    if ID in cExcel.IDs_IP and len(IDs_IP) < n:
                        IDs_IP.append(ID)
                    if ID in cExcel.IDs_family and len(IDs_family) < n:
                        IDs_family.append(ID)
                else:
                    if ID in cExcel.IDs_GO and len(IDs_GO) < n and ID in conditions_dict[condition]:
                        IDs_GO.append(ID)
                    if ID in cExcel.IDs_KW and len(IDs_KW) < n and ID in conditions_dict[condition]:
                        IDs_KW.append(ID)
                    if ID in cExcel.IDs_IP and len(IDs_IP) < n and ID in conditions_dict[condition]:
                        IDs_IP.append(ID)
                    if ID in cExcel.IDs_family and len(IDs_family) < n and ID in conditions_dict[condition]:
                        IDs_family.append(ID)

            with open("%s_GO_leastCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_GO))
            with open("%s_KW_leastCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_KW))
            with open("%s_IP_leastCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_IP))
            with open("%s_family_leastCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs_family))
            IDs_GO_global += IDs_GO
            IDs_KW_global += IDs_KW
            IDs_IP_global += IDs_IP
            IDs_family_global += IDs_family
            intersect_background(IDs_GO_global,"%s_GO_all.txt" %reference_label,"%s_GO_all.txt" %label)
            intersect_background(IDs_KW_global,"%s_KW_all.txt" %reference_label,"%s_KW_all.txt" %label)
            intersect_background(IDs_IP_global,"%s_IP_all.txt" %reference_label,"%s_IP_all.txt" %label)
            intersect_background(IDs_family_global,"%s_family_all.txt" %reference_label,"%s_family_all.txt" %label)
            enrichment_command("%s_GO_mostCV.txt" %label, "%s_GO_all.txt" %label, pwd, "GO")
            enrichment_command("%s_KW_mostCV.txt" %label, "%s_KW_all.txt" %label, pwd, "KW")
            enrichment_command("%s_IP_mostCV.txt" %label, "%s_IP_all.txt" %label, pwd, "IP")
            enrichment_command("%s_family_mostCV.txt" %label, "%s_family_all.txt" %label, pwd, "family")
            enrichment_command("%s_GO_leastCV.txt" %label, "%s_GO_all.txt" %label, pwd, "GO")
            enrichment_command("%s_KW_leastCV.txt" %label, "%s_KW_all.txt" %label, pwd, "KW")
            enrichment_command("%s_IP_leastCV.txt" %label, "%s_IP_all.txt" %label, pwd, "IP")
            enrichment_command("%s_family_leastCV.txt" %label, "%s_family_all.txt" %label, pwd, "family")


    def getSpecificEnrichmentInputAnnotated(self, output_label, n, reference, completeExcel, pwd, target, final_output_label, expressionfiltern=False, annotation_required=False):
        cExcel = CompleteExcel(completeExcel)
        files_mostCV = []
        files_leastCV = []
        if expressionfiltern != False:
            conditions_dict = self.expressionfilter(expressionfiltern)
        for condition in self.conditions:
            IDs = []
            print(len(conditions_dict[condition]))
            label = "%s_%s" %(output_label,"_".join(condition.split()))
            df = self.df.sort_values(condition, ascending=False)
            for j in range(0,df.shape[0]):
                ID = df.iloc[j]["ID"].split(".")[0]
                if expressionfiltern == False:
                    if annotation_required == False and len(IDs) < n:
                        IDs.append(ID)
                    elif annotation_required == True and ID in cExcel.IDs_GO and len(IDs) < n:
                        IDs.append(ID)
                else:
                    if annotation_required == False and len(IDs) < n and ID in conditions_dict[condition]:
                        IDs.append(ID)
                    elif annotation_required == True and ID in cExcel.IDs_GO and len(IDs) < n and ID in conditions_dict[condition]:
                        IDs.append(ID)

            with open("%s_mostCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs))

            IDs_global = IDs
            IDs=[]
            df = self.df.sort_values(condition)
            for j in range(0,df.shape[0]):
                ID = df.iloc[j]["ID"].split(".")[0]
                if expressionfiltern == False:
                    if annotation_required == False and len(IDs) < n:
                        IDs.append(ID)
                    elif annotation_required == True and ID in cExcel.IDs_GO and len(IDs) < n:
                        IDs.append(ID)
                else:
                    if annotation_required == False and len(IDs) < n and ID in conditions_dict[condition]:
                        IDs.append(ID)
                    elif annotation_required == True and ID in cExcel.IDs_GO and len(IDs) < n and ID in conditions_dict[condition]:
                        IDs.append(ID)
            IDs_global += IDs
            with open("%s_leastCV.txt" %label,"wt") as out:
                out.write("\n".join(IDs))

            intersect_background(IDs_global,reference,"%s_all.txt" %label)
            enrichment_command("%s_mostCV.txt" %label, "%s_all.txt" %label, pwd, target)
            enrichment_command("%s_leastCV.txt" %label, "%s_all.txt" %label, pwd, target)

            outname_mostCV = "%s/%s_mostCV_goTerm.tsv" %(target,label.split("/")[-1])
            outname_leastCV = "%s/%s_leastCV_goTerm.tsv" %(target,label.split("/")[-1])
            files_mostCV.append(outname_mostCV)
            files_leastCV.append(outname_leastCV)
        with open ("%s/mostCV_files.txt" %target,"wt") as out:
            for file in files_mostCV:
                out.write("%s\n" %file)
        with open ("%s/leastCV_files.txt" %target,"wt") as out:
            for file in files_leastCV:
                out.write("%s\n" %file)
        enrichment_output_analysis("%s/mostCV_files.txt" %target, final_output_label + "mostCV")
        enrichment_output_analysis("%s/leastCV_files.txt" %target, final_output_label + "leastCV")

class SuperExactTest:
    def __init__(self, file):
        with open(file) as sp:
            self.lines = list(sp)
        self.header = self.lines[0]
        self.name = file.split("/")[-1]

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n","").split(",")
            self.intersections = self.row[0].replace('"','')
            self.intersection_conditions = self.intersections.split(" & ")
            self.degree = self.row[1]
            self.observed_overlap = self.row[2]
            self.expected_overlap = self.row[3]
            self.FE = self.row[4]
            self.Pvalue = self.row[5]
            self.elements = line.split('"')[3].split(", ")
            self.elements = [x for x in self.elements if x]

    def filterByDegree(self,degree_list):
        filtered_lines = [self.header]
        for n in range(1, len(self.lines)):
            line = SuperExactTest.Line(self.lines[n])
            if line.degree in degree_list:
                filtered_lines.append(self.lines[n])
        self.lines = filtered_lines

    def elements_union(self):
        union_list = []
        for n in range(1,len(self.lines)):
            line = SuperExactTest.Line(self.lines[n])
            if line.elements != [""]:
                union_list = union(union_list,line.elements)
        return union_list

    def overlap_stats(self):
        overlap_list = []
        for n in range(1, len(self.lines)):
            line = SuperExactTest.Line(self.lines[n])
            overlap_list.append(int(line.observed_overlap))
        mean_overlap = np.mean(overlap_list)
        median_overlap = np.median(overlap_list)
        return mean_overlap, median_overlap

    def intersection(self, file):
        file = SuperExactTest(file)
        genes1 = self.elements_union()
        genes2 = file.elements_union()
        intersection_genes = intersection(genes1,genes2)
        return intersection_genes

    def stability_analysis(self, output, degree_list=False, file=False,):
        if degree_list:
            self.filterByDegree(degree_list)
        if file:
            sup_file = SuperExactTest(file)
            if degree_list:
                sup_file.filterByDegree(degree_list)
            union_list_1 = self.elements_union()
            mean_overlap_1, median_overlap_1 = self.overlap_stats()
            union_list_2 = sup_file.elements_union()
            mean_overlap_2, median_overlap_2 = sup_file.overlap_stats()
            intersection_list = intersection(union_list_1, union_list_2)
            with open(output, "wt") as out:
                out.write("File\tMean overlap\tMedian overlap\tUnion number\tUnion genes\tFile\t"
                          "Mean overlap\tMedian overlap\tUnion number\tUnion genes\tIntersection number\t"
                          "Intersection genes\n")
                out.write("%s\t%f\t%f\t%i\t%s\t%s\t%f\t%f\t%i\t%s\t%i\t%s\n" %(self.name, mean_overlap_1,
                            median_overlap_1, len(union_list_1), ",".join(union_list_1), sup_file.name,
                            mean_overlap_2, median_overlap_2, len(union_list_2), ",".join(union_list_2),
                            len(intersection_list), ",".join(intersection_list)))
        else:
            union_list = self.elements_union()
            mean_overlap, median_overlap = self.overlap_stats()
            with open(output, "wt") as out:
                out.write("Mean overlap\tMedian overlap\tUnion number\tUnion genes\n")
                out.write("%f\t%f\t%i\t%s\n" %(mean_overlap, median_overlap, len(union_list), ",".join(union_list)))

    def initial_sizes(self):
        initial_size_dict = {}
        for k in range(1,len(self.lines)):
            line = self.lines[k].split(",")
            if line[1] == "1":
                condition = line[0].replace('"','')
                initial_size = line[2]
                initial_size_dict[condition] = int(initial_size)
        return initial_size_dict

    def DE_analysis(self, output):
        times = ["12h/24h", "24h/48h", "48h/72h", "72h/96h", "96h"]
        sametime_lines = [self.header]
        timelapse_lines = [self.header]
        initial_size_dict = self.initial_sizes()
        hour = []
        lapse = []
        hour_2 = []
        exposure = []


        self.filterByDegree(["2"])
        for k in range(1, len(self.lines)):
            line = SuperExactTest.Line(self.lines[k])
            line.times = []
            line.exposures = []
            for condition in line.intersection_conditions:
                line.times.append(condition.split("_")[2].replace('"',''))
                if condition.split("_")[0].split("de")[1] == "4":
                    line.exposures.append("first")
                elif condition.split("_")[0].split("de")[1] == "5":
                    line.exposures.append("second")
            ### Gene program delay analysis
            if int(line.times[0]) < 5 and line.times[0] == line.times[1] and line.exposures[0] != line.exposures[1]:
                sametime_lines.append(line)
                for n in range(0,int(line.observed_overlap)):
                    hour.append(times[int(line.times[0])-1])
                    lapse.append("Not delayed")
            if int(line.times[0]) < 5 and int(line.times[0]) + 1 == int(line.times[1]) and line.exposures[0] != line.exposures[1]:
                timelapse_lines.append(line)
                for n in range(0, int(line.observed_overlap)):
                    hour.append(times[int(line.times[0])-1])
                    lapse.append("Second exposure delayed")
            if int(line.times[1]) < 5 and int(line.times[1]) + 1 == int(line.times[0]) and line.exposures[0] != line.exposures[1]:
                timelapse_lines.append(line)
                for n in range(0, int(line.observed_overlap)):
                    hour.append(times[int(line.times[1])-1])
                    lapse.append("First exposure delayed")
            ### Timeline analysis
            if int(line.times[0]) < 5 and int(line.times[0]) + 1 == int(line.times[1]) and line.exposures[0] == line.exposures[1]:
                timelapse_lines.append(line)
                if line.exposures[0] == "first":
                    for n in range(0, int(line.observed_overlap)):
                        hour_2.append(times[int(line.times[0])-1])
                        exposure.append("First exposure")
                elif line.exposures[0] == "second":
                    for n in range(0, int(line.observed_overlap)):
                        hour_2.append(times[int(line.times[0])-1])
                        exposure.append("Second exposure")

        hour.reverse()
        lapse.reverse()
        df_dict = dict()
        df_dict["Overlap"] = hour
        df_dict["Delay"] = lapse

        df = pd.DataFrame.from_dict(df_dict)

        fig = plt.figure()
        sns.set_theme()
        sns.countplot(x="Overlap", hue="Delay",  data=df)
        fig.savefig(output)
        plt.close()

        hour_2.reverse()
        exposure.reverse()
        df_dict = dict()
        df_dict["Overlap"] = hour_2
        df_dict["Exposure"] = exposure

        df = pd.DataFrame.from_dict(df_dict)
        fig = plt.figure()
        sns.set_theme()
        sns.countplot(x="Overlap", hue="Exposure", data=df)
        fig.savefig("%s_timeline.png" %output.split(".")[0])

        times = ["12h/24h", "24h/48h", "48h/72h", "72h/96h", "96h"]
        sametime_lines = [self.header]
        timelapse_lines = [self.header]
        hour = []
        lapse = []
        hour_2 = []
        exposure = []

    ##### Normalized by the highest overlaping it can have.
        for k in range(1, len(self.lines)):
            line = SuperExactTest.Line(self.lines[k])
            line.times = []
            line.exposures = []
            for condition in line.intersection_conditions:
                line.times.append(condition.split("_")[2].replace('"',''))
                if condition.split("_")[0].split("de")[1] == "4":
                    line.exposures.append("first")
                elif condition.split("_")[0].split("de")[1] == "5":
                    line.exposures.append("second")
            max_overlap = min([initial_size_dict[line.intersection_conditions[0]],
                               initial_size_dict[line.intersection_conditions[1]]])
            normalized_overlap = round(int(line.observed_overlap) / max_overlap * 100)
            ### Gene program delay analysis
            if int(line.times[0]) < 5 and line.times[0] == line.times[1] and line.exposures[0] != line.exposures[1]:
                sametime_lines.append(line)
                for n in range(0,normalized_overlap):
                    hour.append(times[int(line.times[0])-1])
                    lapse.append("Not delayed")
            if int(line.times[0]) < 5 and int(line.times[0]) + 1 == int(line.times[1]) and line.exposures[0] != line.exposures[1]:
                timelapse_lines.append(line)
                for n in range(0, normalized_overlap):
                    hour.append(times[int(line.times[0])-1])
                    lapse.append("Second exposure delayed")
            if int(line.times[1]) < 5 and int(line.times[1]) + 1 == int(line.times[0]) and line.exposures[0] != line.exposures[1]:
                timelapse_lines.append(line)
                for n in range(0, normalized_overlap):
                    hour.append(times[int(line.times[1])-1])
                    lapse.append("First exposure delayed")
            ### Timeline analysis
            if int(line.times[0]) < 5 and int(line.times[0]) + 1 == int(line.times[1]) and line.exposures[0] == line.exposures[1]:
                timelapse_lines.append(line)
                if line.exposures[0] == "first":
                    for n in range(0, normalized_overlap):
                        hour_2.append(times[int(line.times[0])-1])
                        exposure.append("First exposure")
                elif line.exposures[0] == "second":
                    for n in range(0, normalized_overlap):
                        hour_2.append(times[int(line.times[0])-1])
                        exposure.append("Second exposure")

        hour.reverse()
        lapse.reverse()
        df_dict = dict()
        df_dict["Overlap"] = hour
        df_dict["Delay"] = lapse

        df = pd.DataFrame.from_dict(df_dict)

        fig = plt.figure()
        sns.set_theme()
        sns.countplot(x="Overlap", hue="Delay",  data=df)
        plt.ylabel("Percentage with respect to maximum overlap")
        fig.savefig("%s_percentage.png" %output.split(".")[0])
        plt.close()

        hour_2.reverse()
        exposure.reverse()
        df_dict = dict()
        df_dict["Overlap"] = hour_2
        df_dict["Exposure"] = exposure

        df = pd.DataFrame.from_dict(df_dict)
        fig = plt.figure()
        sns.set_theme()
        sns.countplot(x="Overlap", hue="Exposure", data=df)
        plt.ylabel("Percentage with respect to maximum overlap")
        fig.savefig("%s_timeline_percentage.png" %output.split(".")[0])

    def DE_special_cases_MG(self, genes_file, output_label):
        # getting genes dictionary
        genes_dict = dict()
        conditions = []
        with open(genes_file) as gf:
            gf = list(gf)

        for condition in gf[0].split("\t"):
            genes_dict[condition] = []
            conditions.append(condition)


        for j in range(1, len(gf)):
            for k in range(0, len(gf[j].split("\t"))):
                genes_dict[conditions[k]].append(gf[j].split("\t")[k])

        self.filterByDegree(["2"])
        ### getting overlap dict
        overlap_dict = dict()
        for j in range(1, len(self.lines)):
            line = SuperExactTest.Line(self.lines[j])
            overlap_dict[line.intersections] = line.elements
        ### Overlap between DE not overlaping in first exposure 12hvs24h and and second exposure 24h vs 48h

        de4_12vs24 = union(genes_dict["de4_MG_1"], genes_dict["de4_MG_2"])
        de4_notoverlapping_12vs24 = []
        for gene in de4_12vs24:
            if gene not in overlap_dict["de4_MG_1 & de4_MG_2"]:
                de4_notoverlapping_12vs24.append(gene)

        de5_24vs48 = union(genes_dict["de5_MG_2"], genes_dict["de5_MG_3"])
        de5_notoverlapping_24vs48 = []
        for gene in de5_24vs48:
            if gene not in overlap_dict["de5_MG_2 & de5_MG_3"]:
                de5_notoverlapping_24vs48.append(gene)

        not_overlapping = union(de4_notoverlapping_12vs24, de5_notoverlapping_24vs48)

        result = intersection(de4_notoverlapping_12vs24, de5_notoverlapping_24vs48)

        with open("%s_transcripts.txt" % output_label, "wt") as out:
            for gene in result:
                out.write("%s\n" %gene)

        with open("%s_stats.txt" %output_label, "wt") as out:
            out.write("genes in union de4_1 and de4_2\tgenes in union de5_2 and de5_3\tgenes not overlapping de4_1 vs"
                      " de4_2\tgenes not overlapping de5_2 vs de5_3\tgenes in union\tgenes in intersection\tratio\n")
            out.write("%i\t%i\t%i\t%i\t%i\t%i\t%f\n" %(len(de4_12vs24), len(de5_24vs48), len(de4_notoverlapping_12vs24),
                                                       len(de5_notoverlapping_24vs48), len(not_overlapping), len(result),
                                                       len(result)/len(not_overlapping)))

    def DE_special_cases_SG(self, genes_file, output_label):
        # getting genes dictionary
        genes_dict = dict()
        conditions = []
        with open(genes_file) as gf:
            gf = list(gf)

        for condition in gf[0].split("\t"):
            genes_dict[condition] = []
            conditions.append(condition)

        for j in range(1, len(gf)):
            for k in range(0, len(gf[j].split("\t"))):
                genes_dict[conditions[k]].append(gf[j].split("\t")[k])

        self.filterByDegree(["2"])
        ### getting overlap dict
        overlap_dict = dict()
        for j in range(1, len(self.lines)):
            line = SuperExactTest.Line(self.lines[j])
            overlap_dict[line.intersections] = line.elements
        
        ### Overlap between DE that appears as overlapping for every comparison - First exposure
        intersection_list = intersection(overlap_dict["de4_SG_1 & de4_SG_2"], overlap_dict["de4_SG_2 & de4_SG_3"])
        intersection_list = intersection(intersection_list, overlap_dict["de4_SG_3 & de4_SG_4"])
        intersection_list = intersection(intersection_list, overlap_dict["de4_SG_4 & de4_SG_5"])

        union_list = union(overlap_dict["de4_SG_1 & de4_SG_2"], overlap_dict["de4_SG_2 & de4_SG_3"])
        union_list = union(union_list, overlap_dict["de4_SG_3 & de4_SG_4"])
        union_list = union(union_list, overlap_dict["de4_SG_4 & de4_SG_5"])

        with open("%s_SG_SpecialCase_OverlappingForEveryComparison_FirstExposure_transcripts.txt" %output_label, "wt") as out:
            for gene in intersection_list:
                out.write("%s\n" %gene)

        with open("%s_SG_SpecialCase_OverlappingForEveryComparison_FirstExposure_stats.txt" %output_label, "wt") as out:
            out.write("union of genes for every comparison\tgenes overlapping in every comparison\tratio\n")
            out.write("%i\t%i\t%f\n" %(len(union_list), len(intersection_list), len(intersection_list)/len(union_list)))

        ### Overlap between DE that appears as overlapping for every comparison - Second exposure
        intersection_list = intersection(overlap_dict["de5_SG_1 & de5_SG_2"], overlap_dict["de5_SG_2 & de5_SG_3"])
        intersection_list = intersection(intersection_list, overlap_dict["de5_SG_3 & de5_SG_4"])
        intersection_list = intersection(intersection_list, overlap_dict["de5_SG_4 & de5_SG_5"])
        intersection_list_second = intersection_list

        union_list = union(overlap_dict["de5_SG_1 & de5_SG_2"], overlap_dict["de5_SG_2 & de5_SG_3"])
        union_list = union(union_list, overlap_dict["de5_SG_3 & de5_SG_4"])
        union_list = union(union_list, overlap_dict["de5_SG_4 & de5_SG_5"])

        with open("%s_SG_SpecialCase_OverlappingForEveryComparison_SecondExposure_transcripts.txt" %output_label, "wt") as out:
            for gene in intersection_list:
                out.write("%s\n" % gene)

        with open("%s_SG_SpecialCase_OverlappingForEveryComparison_SecondExposure_stats.txt" %output_label, "wt") as out:
            out.write("union of genes for every comparison\tgenes overlapping in every comparison\tratio\n")
            out.write(
                "%i\t%i\t%f\n" % (len(union_list), len(intersection_list), len(intersection_list) / len(union_list)))

        ###DE that don't overlap in at least one comparison of First and that don't overlap in every comparison of second

        de4_12vs24 = union(genes_dict["de4_SG_1"], genes_dict["de4_SG_2"])
        de4_24vs48 = union(genes_dict["de4_SG_2"], genes_dict["de4_SG_3"])
        de4_48vs72 = union(genes_dict["de4_SG_3"], genes_dict["de4_SG_4"])
        de4_72vs96 = union(genes_dict["de4_SG_4"], genes_dict["de4_SG_5"])

        de4_notoverlapping_12vs24 = intersection(de4_12vs24, overlap_dict["de4_SG_1 & de4_SG_2"], reverse=True)
        de4_notoverlapping_24vs48 = intersection(de4_24vs48, overlap_dict["de4_SG_2 & de4_SG_3"], reverse=True)
        de4_notoverlapping_48vs72 = intersection(de4_48vs72, overlap_dict["de4_SG_3 & de4_SG_4"], reverse=True)
        de4_notoverlapping_72vs96 = intersection(de4_72vs96, overlap_dict["de4_SG_4 & de4_SG_5"], reverse=True)

        union_notoverlapping = union(de4_notoverlapping_12vs24, de4_notoverlapping_24vs48)
        union_notoverlapping = union(union_notoverlapping, de4_notoverlapping_48vs72)
        union_notoverlapping = union(union_notoverlapping, de4_notoverlapping_72vs96)

        result = intersection(union_notoverlapping, intersection_list_second, reverse=True)

        with open("%s_SG_SpecialCase_NotOverlappingInFirstNotOverlappingEveryComparisonSecond_transcripts.txt" %output_label, "wt") as out:
            for gene in result:
                out.write("%s\n" %gene)

        with open("%s_SG_SpecialCase_NotOverlappingInFirstNotOverlappingEveryComparisonSecond_stats.txt" %output_label, "wt") as out:
            out.write("union of genes\tgenes that match the case\tratio\n")
            out.write("%i\t%i\t%f\n" %(len(union_notoverlapping), len(result), len(result)/len(union_notoverlapping)))

class MaSigProClusters:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(fl)
        self.dict = self.clusters2dict()

    def clusters2dict(self):
        cluster_dict = dict()
        for k in range(1, len(self.lines)):
            row = self.lines[k].replace("\n", "").split("\t")
            name = row[0].replace('"', '')
            cluster = row[1]
            if cluster not in cluster_dict:
                cluster_dict[cluster] = [name]
            else:
                cluster_dict[cluster].append(name)
        return cluster_dict

    def get_gene_lists(self, out_label):
        for key in self.dict:
            with open("%s_cluster%s.txt" %(out_label, key), "wt") as out:
                out.write("\n".join(self.dict[key]))

    def GO_annotation(self, background, out_label):
        for key in self.dict:
            out_name = "%s_cluster%s.txt" %(out_label, key)
            with open(out_name, "wt") as out:
                out.write("\n".join(self.dict[key]))
            os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s" % (
            out_name, background))

    def clusters2table(self):
        x="FUTURE DEVELOPMENT"

class GOBackground:
    def __init__(self,file):
        with open(file) as fl:
            self.lines = list(fl)
        self.dict = self.get_dict()

    def get_dict(self):
        go_dict = dict()
        for line in self.lines:
            id=line.split("\t")[0]
            go = line.replace("\n","").split("\t")[1]
            if "," in go:
                go = go.split(",")
            else:
                go = [go]
            go_dict[id] = go
        return go_dict

    def get_ancestors(self,obo,output):
        godag = GODag(obo, optional_attrs={'relationship'})
        optional_relationships = {'part_of','regulates', 'negatively_regulates', 'positively_regulates'}
        b_dict = dict()
        for key in self.dict:
            b_dict[key] = []
            for go in self.dict[key]:
                b_dict[key].append(go)
                gosubdag_r0 = GoSubDag([go], godag, relationships=optional_relationships, prt=None)
                if go in gosubdag_r0.rcntobj.go2ancestors:
                    b_dict[key] = b_dict[key] + list(gosubdag_r0.rcntobj.go2ancestors[go])
            b_dict[key] = deduplicate_list(b_dict[key])
        with open(output,"wt") as out:
            for key in b_dict:
                out.write("%s\t%s\n" %(key, ",".join(b_dict[key])))

class StringDBGO:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(fl)
        self.fdr_dict, self.strength_dict = self.get_dict()

    class Line:
        def __init__(self,line):
            self.row = line.replace("\n","").split("\t")
            self.id = self.row[0]
            self.description = self.row[1]
            self.fdr = self.row[5]
            self.strength = self.row[4]


    def get_dict(self):
        fdr_dict = dict()
        strength_dict = dict()
        for k in range(1,len(self.lines)):
            ln = StringDBGO.Line(self.lines[k])
            fdr_dict[ln.description] = -1* math.log(float(ln.fdr), 10)
            strength_dict[ln.description] = float(ln.strength)
        return fdr_dict,strength_dict

class Bed:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(file)

    class Line:
        def __init__(self,line):
            self.row = line.replace("\n",""). split("\t")
            self.chr = self.row[0]
            self.init = self.row[1]
            self.end = self.row[2]
            self.id = self.row[3]
            self.strand = self.row[5]
            self.flag = self.row[6]
            self.cigar = self.row[7]

    def intersect(self, i_list, output):
        with open(output) as out:
            for line in self.lines:
                bed_line = Bed.Line(line)
                if bed_line.id in i_list:
                    out.write(line)

class IxoriDBOut:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(fl)
        self.ids = self.get_ids()
        self.protein_dict = self.get_proteins()
        self.cds_dict = self.get_CDS()

    class Line:
        def __init__(self, line):
            self.row = line.replace('"','').replace("\n","").split(",")
            self.id = self.row[0]
            self.protein = self.row[-1]
            self.cds = self.row[-2]

    def get_ids(self):
        ids = []
        for k in range(1, len(self.lines)):
            line = IxoriDBOut.Line(self.lines[k])
            ids.append(line.id)
        return ids

    def get_proteins(self):
        protein_dict = dict()
        for k in range(1, len(self.lines)):
            line = IxoriDBOut.Line(self.lines[k])
            protein_dict[line.id] = line.protein
        return protein_dict

    def get_CDS(self):
        CDS_dict = dict()
        for k in range(1, len(self.lines)):
            line = IxoriDBOut.Line(self.lines[k])
            CDS_dict[line.id] = line.cds
        return CDS_dict

class Outfmt6:
    def __init__(self,file):
        with open(file) as fl:
            self.lines = list(fl)
        self.pname_dict = self.get_pname_dict()

    class Line:
        def __init__(self,line):
            self.row = line.replace("\n","").split("\t")
            self.id = self.row[0]
            self.pname = self.row[1]
            if "|" in self.pname:
                self.pname = self.pname.split("|")[1]
            self.bitscore = float(self.row[-1].replace("\n",""))

    def getnames(self):
        names = []
        for n in range(0, len(self.lines)):
            line = Outfmt6.Line(self.lines[n])
            names.append(line.pname)
        names = deduplicate_list(names)
        return names

    def get_pname_dict(self):
        name_dict = dict()
        for line in self.lines:
            line = Outfmt6.Line(line)
            if line.id not in name_dict:
                name_dict[line.id] = line.pname
        return name_dict

    def get_pname_bitscore_dict(self, bitscore=0):
        name_dict = dict()
        for line in self.lines:
            line = Outfmt6.Line(line)
            if line.id not in name_dict and line.bitscore > bitscore:
                name_dict[line.id] = [line.pname]
            elif line.id in name_dict and line.bitscore > bitscore:
                name_dict[line.id].append(line.pname)
        return name_dict

    def family_annotation(self, annotation_file, output):
        with open(annotation_file) as af:
            af = list(af)
        family_dict = dict()
        for line in af:
            row = line.replace("\n","").split("\t")
            pname = row[0]
            family = row[8]
            family_dict[pname] = family
        with open(output, "wt") as out:
            for id in self.pname_dict:
                if self.pname_dict[id] in family_dict:
                    out.write("%s\t%s\t%s\n" %(id, self.pname_dict[id], family_dict[self.pname_dict[id]]))

    def get_sequences(self, fasta, output, bitscore=0):
        id_dict = self.get_pname_bitscore_dict(bitscore)
        fa = Fasta(fasta)

        ## list
        id_list = []

        for key in id_dict:
            if isinstance(id_dict[key], list):
                id_list = id_list + id_dict[key]
            else:
                id_list.append(id_dict[key])


        fa.filter(id_list, output)

class BlatGff:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(fl)
            self.coord_dict = self.get_coords()

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.chr = self.row[0]
            self.program = self.row[1]
            self.type = self.row[2]
            self.init = self.row[3]
            self.end = self.row[4]
            self.score = self.row[5]
            self.strand = self.row[6]
            self.parent = self.row[8].split("Parent=")[1].split(";")[0]
            self.target = self.row[8].split("Target=")[1].split(";")[0]

    def get_coords(self):
        coord_dict = dict()
        for k in range(1, len(self.lines)):
            line = BlatGff.Line(self.lines[k])
            if line.strand == "++":
                if line.parent in coord_dict:
                    coord_dict[line.parent][1] = line.end
                else:
                    coord_dict[line.parent] = [line.init, line.end]
            else:
                if line.parent in coord_dict:
                    coord_dict[line.parent][1] = line.end
                else:
                    coord_dict[line.parent] = [line.init, line.end]

        return coord_dict

    def depurate_gff(self, reference, output):

        parent_list = []

        ref = Fasta(reference)
        with open(output,"wt") as out:
            out.write(self.lines[0])
            for k in range(1,len(self.lines)):
                line = BlatGff.Line(self.lines[k])

                if line.strand == "++":
                    line.strand = "+"
                    parent_start = int(self.coord_dict[line.parent][0])
                    parent_end = int(self.coord_dict[line.parent][1])
                    coord_start = int(line.init)
                    coord_end = int(line.end)
                elif line.strand == "+-":
                    line.strand = "-"
                    parent_end = ref.len_dict[line.chr] - int(self.coord_dict[line.parent][0])
                    parent_start = ref.len_dict[line.chr] - int(self.coord_dict[line.parent][1])
                    coord_end = ref.len_dict[line.chr] - int(line.init)
                    coord_start = ref.len_dict[line.chr] - int(line.end)

                if line.parent not in parent_list:
                    out.write("%s\t%s\tgene\t%i\t%i\t%s\t%s\t."
                              "\tID=%s;gene_id=%s;transcript_id=%s.1\n" %(line.chr, line.program, parent_start,
                                                                        parent_end, line.score, line.strand,
                                                                        line.parent, line.parent, line.parent))
                    out.write("%s\t%s\ttranscript\t%i\t%i\t%s\t%s\t."
                              "\tID=%s.1;Parent=%s;transcript_id=%s.1\n" % (
                              line.chr, line.program, parent_start,
                              parent_end, line.score, line.strand,
                              line.parent, line.parent, line.parent))
                    out.write("%s\t%s\texon\t%i\t%i\t%s\t%s\t."
                              "\tParent=%s;Target=%s\n" % (
                                  line.chr, line.program, coord_start,
                                  coord_end, line.score, line.strand,
                                  line.parent, line.target))
                    parent_list.append(line.parent)
                else:
                    out.write("%s\t%s\texon\t%i\t%i\t%s\t%s\t."
                              "\tParent=%s;Target=%s\n" % (
                                  line.chr, line.program, coord_start,
                                  coord_end, line.score, line.strand,
                                  line.parent, line.target))

class GffCompareTrackingGroup:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(fl)
        self.transcripts = self.get_transcripts()

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.transcript = line.split("q1:")[1].split("|")[0]

    def get_transcripts(self):
        transcripts = []
        for k in range(1, len(self.lines)):
            line = GffCompareTrackingGroup.Line(self.lines[k])
            transcripts.append(line.transcript)
        return transcripts

class Gff:
    def __init__(self, file):
        with open(file) as fl:
            self.lines = list(fl)
        self.gene_coords = self.get_gene_coords()

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.chr = self.row[0]
            self.program = self.row[1]
            self.type = self.row[2]
            self.init = self.row[3]
            self.end = self.row[4]
            self.score = self.row[5]
            self.strand = self.row[6]
            self.frame = self.row[7]
            self.id = self.row[8].split("ID=")[1].split(";")[0]
            self.description_noID = self.row[8].replace("ID=%s;" %self.id, "")
            self.description = self.row[8]

    def get_ids(self, label):
        ids = []
        for line in self.lines:
            if not line.startswith("#"):
                line = Gff.Line(line)
                if line.type == label:
                    ids.append(line.id)
        return ids

    def get_sequences_from_IDs(self, label, fasta, output):
        fa = Fasta(fasta)
        ids = self.get_ids(label)
        fa.filter(ids, output)

    def get_gene_coords(self):
        coord_dict = dict()
        for line in self.lines:
            if not line.startswith("#"):
                if "gene" in line:
                    ln = Gff.Line(line)
                    coord_dict[ln.id] = [ln.chr, ln.strand, ln.init, ln.end]

    def filter_by_score(self, score, output):
        with open(output, "wt") as out:
            for line in self.lines:
                if not line.startswith("#"):
                    ln = Gff.Line(line)
                    if float(ln.score) > score:
                        out.write(line)

    def add_CDS(self, output):
        with open(output, "wt") as out:
            exon = []
            for line in self.lines:
                if "\texon\t" in line:
                    exon.append(line)
                else:
                    if exon != []:
                        out.write("".join(exon))
                        out.write("".join(exon).replace("exon", "CDS"))
                        exon = []
                    out.write(line)

    def add_exon(self,output):
        with open(output, "wt") as out:
            for line in self.lines:
                ln = Gff.Line(line)
                out.write(
                    "%s\t%s\tgene\t%s\t%s\t%s\t%s\t.\tID=%s;%s\n" % (ln.chr, ln.program, ln.init, ln.end,
                                                                                ln.score, ln.strand, ln.id, ln.description_noID))
                out.write(
                    "%s\t%s\ttranscript\t%s\t%s\t%s\t%s\t.\tID=%s.1;Parent=%s;%s\n" % (ln.chr, ln.program, ln.init, ln.end,
                                                                                ln.score, ln.strand, ln.id, ln.id, ln.description_noID))
                out.write("%s\t%s\texon\t%s\t%s\t%s\t%s\t.\tID=%s.1.exon1;Parent=%s.1;%s\n" %(ln.chr, ln.program, ln.init, ln.end,
                                                                                     ln.score, ln.strand, ln.id, ln.id, ln.description_noID))

    def map_to_fasta(self, fasta, fasta_ref, out_label):
        fa = Fasta(fasta)
        fa_ref = Fasta(fasta_ref)
        write = False

        with open("%s_mapped.gff" %out_label, "wt") as out_mapped:
            with open("%s_unmapped.fa" %out_label, "wt") as out_unmapped:
                for line in self.lines:
                    if line.startswith("#"):
                        out_mapped.write(line)
                    else:
                        gff_line = Gff.Line(line)
                        if gff_line.type == "gene":
                            write=False
                            id = gff_line.id + "-PA"
                            if id in fa.dict.keys():
                                write=True
                                out_mapped.write(line)
                            else:
                                out_unmapped.write(">%s\n%s\n" %(id, fa_ref.dict[id]))
                        else:
                            if write == True:
                                out_mapped.write(line)

    def to_tsv(self, type, output):
        with open(output, "wt") as out:
            out.write("ID\tChr\tSoftware\tClass\tStart\tEnd\tScore\tStrand\tFrame\tDescription\n")
            for line in self.lines:
                if not line.startswith("#"):
                    gff_line = Gff.Line(line)
                    if gff_line.type == type:
                        out.write("%s\t%s" %(gff_line.id))

class SampleSheet:
    def __init__(self,file):
        with open(file) as fl:
            self.lines = list(fl)
        self.dict = self.get_dict()

    def get_dict(self):
        group_dict = dict()
        for line in self.lines:
            row = getrow(line)
            group_dict[row[2]] = row[1]
        return group_dict



def correct_rfam_gff(gff_file, output):
    gff = Gff(gff_file)
    with open(output, "wt") as out:
        for line in gff.lines:
            ln = Gff.Line(line)
            if ln.strand == "-":
                if ln.program == "infernal":
                    start = ln.end
                    end = ln.init
                else:
                    start = ln.init
                    end = ln.end
            else:
                start = ln.init
                end = ln.end
            if "gene" in line:
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(ln.chr, ln.program, ln.type, start, end, ln.score, ln.strand, ln.frame, ln.description))
            elif "transcript" in line:
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(ln.chr, ln.program, ln.type, start, end, ln.score, ln.strand, ln.frame, ln.description))
                out.write("%s\t%s\texon\t%s\t%s\t%s\t%s\t%s\tID=%s.exon1;Parent=%s;%s\n" % (
                ln.chr, ln.program, start, end, ln.score, ln.strand, ln.frame, ln.id, ln.id, ";".join(ln.description.split("Parent")[1].split(";")[1::])))


def ires_analysis(gff, ref_gff, output):
    ref_gff = Gff(ref_gff)
    gff = Gff(gff)
    with open(output, "wt") as out:
        out.write("Chr\tProgram\tClass\tStart\tEnd\tScore\tStrand\tFrame\tDescription\tGene\tFull name\tCoordinates\tLocation\tMinimum Distance\n")
        for line in gff.lines:
            if "gene" in line and "Iron_response" in line:
                ln = Gff.Line(line)
                coords = [ln.chr, ln.strand, int(ln.init), int(ln.end)]
                closest_gene, fullname, gene_coord, location, distance = find_adjacent_gene(coords, ref_gff)
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(line.replace("\n", ""), closest_gene, fullname, gene_coord, location, distance))



def find_adjacent_gene(coords, ref_gff):
    closest_gene = None
    min_distance = float('inf')
    for line in ref_gff.lines:
        if "gene" in line:
            ln = Gff.Line(line)
            if ln.chr == coords[0] and ln.strand == coords[1]:
                distance = min(abs(coords[2] - int(ln.init)), abs(coords[3] - int(ln.end)))
                if distance < min_distance:
                    min_distance = distance
                    closest_gene = ln.id
                    if "full_name" in line:
                        fullname = line.split("full_name=")[1].split(";")[0]
                    else:
                        fullname = closest_gene
                    gene_coord = "%s:%s" %(ln.init, ln.end)
    for line in ref_gff.lines:
        if closest_gene in line:
            ln = Gff.Line(line)
            if coords[1] == "-" and coords[2] < int(ln.end) and coords[3] > int(ln.init):
                location = ln.type
                min_distance = 0
            elif coords[1] == "+" and coords[2] > int(ln.init) and coords[3] < int(ln.end):
                location = ln.type
                min_distance = 0
            else:
                location= "adjacent"
    if location == "mRNA" or location == "gene" or location == "transcript":
        location = "intronic"

    return closest_gene, fullname, gene_coord, location, min_distance



def crossreference_OGS1(cystatin, kunitz, gff, output):
    ogs = Gff(gff)
    cys = Fasta(cystatin)
    kun = Fasta(kunitz)

    symbols = list(cys.dict.keys()) + list(kun.dict.keys())
    annotated_symbols = []
    with open(output, "wt") as out:
        out.write("ID\tSymbol\tFull name\n")
        for line in ogs.lines:
            if "symbol" in line:
                symbol = line.split("symbol=")[1].replace("\n","")
                if ";" in symbol:
                    symbol = symbol.split(";")[0]
                if symbol in symbols and symbol not in annotated_symbols:
                    full_name = line.split("full_name=")[1].split(";")[0]
                    id = line.split("ID=")[1].split(";")[0]
                    out.write("%s\t%s\t%s\n" %(id, symbol, full_name))
                    annotated_symbols.append(symbol)

def ncRNA_add_type(gff, output):
    gff = Gff(gff)

    with open(output, "wt") as out:
        for line in gff.lines:
            gff_line = Gff.Line(line)
            pos_id = "_%s_%s" %(gff_line.chr.replace("Iricinus_Scaffold_",""), gff_line.init)
            type = gff_line.id.split(pos_id)[0]
            if "mir" in type.lower():
                type = "miRNA"
            if line.replace("\n","")[-1] != ";":
                out.write(line.replace("\n","") + ";type=%s;\n" %type)
            else:
                out.write(line.replace("\n", "") + "type=%s;\n" % type)

def lncRNA_gff_rename(gff, consensus, output):
    gff = Gff(gff)
    consensus = Fasta(consensus)

    ids = consensus.dict.keys()

    with open(output, "wt") as out:
        tname_counter ={}
        for line in gff.lines:
            if not line.startswith("#"):
                gff_line = Gff.Line(line)
                if gff_line.type == "mRNA":
                    consensus_id = gff_line.id.split(".")[0]
                    scaff_n = gff_line.chr.replace("Iricinus_Scaffold_", "")
                    if consensus_id in ids:
                        tname = "HC_lncRNA_%s_%s" %(scaff_n, gff_line.init)
                        if tname in tname_counter:
                            new_tname = tname + ".%i" %tname_counter[tname]
                            tname_counter[tname] = tname_counter[tname] + 1
                            tname = new_tname
                        else:
                            tname_counter[tname] = 2
                    else:
                        tname = "LC_lncRNA_%s_%s" %(scaff_n, gff_line.init)
                        if tname in tname_counter:
                            new_tname = tname + ".%i" %tname_counter[tname]
                            tname_counter[tname] = tname_counter[tname] + 1
                            tname = new_tname
                        else:
                            tname_counter[tname] = 2
                    out.write("%s\t%s\tgene\t%s\t%s\t%s\t%s\t.\tID=%s;Name=%s;\n" %(gff_line.chr, gff_line.program, gff_line.init,
                                                                                     gff_line.end, gff_line.score, gff_line.strand,
                                                                                    tname, tname))
                    out.write("%s\t%s\ttranscript\t%s\t%s\t%s\t%s\t.\tID=%s.1;Parent=%s;\n" %(gff_line.chr, gff_line.program, gff_line.init,
                                                                                     gff_line.end, gff_line.score, gff_line.strand,
                                                                                    tname, tname))
                elif gff_line.type == "exon":
                    exon = gff_line.id.split(".")[-1]
                    out.write("%s\t%s\texon\t%s\t%s\t%s\t%s\t.\tID=%s.1.%s;Parent=%s.1;\n" % (
                    gff_line.chr, gff_line.program, gff_line.init,
                    gff_line.end, gff_line.score, gff_line.strand,
                    tname, exon, tname))
            else:
                out.write(line)

def get_trna_gff(trnascan, output):
    with open(trnascan) as fl:
        lines = list(fl)
    with open(output, "wt") as out:
        for line in lines:
            if not line.startswith("#"):
                row = line.replace("\n","").split("\t")
                scaff = row[0]
                software = row[1]
                type = row[2]
                start = row[3]
                end = row[4]
                score= row[5]
                strand = row[6]
                frame = row[7]
                description = row[8]

                if type == "tRNA":
                    isotype = description.split("isotype=")[1].split(";")[0]
                    anticodon = description.split("anticodon=")[1].split(";")[0]
                    gene_biotype = description.split("gene_biotype=")[1].split(";")[0]
                    ID = description.split("ID=")[1].split(";")[0]

                    out.write("%s\t%s\tgene\t%s\t%s\t%s\t%s\t%s\tID=tRNA_%s_%s;isotype=%s;anticodon=%s;type=%s;\n" %(scaff, software, start, end, score, strand,
                                                                 frame, scaff, start, isotype, anticodon, gene_biotype ))
                    out.write(
                        "%s\t%s\ttranscript\t%s\t%s\t%s\t%s\t%s\tID=tRNA_%s_%s.1;Parent=tRNA_%s_%s;isotype=%s;anticodon=%s;type=%s;\n" % (
                        scaff, software, start, end, score, strand,
                        frame, scaff, start, scaff, start,isotype, anticodon, gene_biotype ))
                    out.write(
                        "%s\t%s\texon\t%s\t%s\t%s\t%s\t%s\tID=tRNA_%s_%s.1.exon1;Parent=tRNA_%s_%s.1;isotype=%s;anticodon=%s;type=%s;\n" % (
                            scaff, software, start, end, score, strand,
                            frame, scaff, start, scaff, start,isotype, anticodon, gene_biotype ))

def excel_cystatin(file, output):
    fa = Fasta(file)
    with open(output,"wt") as out:
        out.write("Protein ID\tProtein description\tScaffold\tStart\tEnd\tStrand\tLength\tSequence\n")
        for key in fa.full_dict:
            id = key.split(" ")[0].replace(">","")
            description = " ".join(key.split(" ")[1:])
            scaffold = key.split("[")[1].split(":")[0]
            start = key.split("[")[1].split(":")[1].split("-")[0]
            end = key.split("[")[1].split(":")[1].split("-")[1].split(" ")[0]
            strand = key.split("[")[1].split(" ")[1]
            length = len(fa.full_dict[key])
            seq = fa.full_dict[key]
            out.write("%s\t%s\t%s\t%s\t%s\t'%s\t%s\t%s\n" %(id, description, scaffold, start, end,
                                                               strand, length, seq))

def excel_kunitz(file, output):
    fa = Fasta(file)
    with open(output,"wt") as out:
        out.write("Protein ID\tProtein description\tScaffold\tStart\tEnd\tStrand\tKunitz domains\tLength\tSequence\n")
        for key in fa.full_dict:
            id = key.split(" ")[0].replace(">","")
            description = " ".join(key.split(" ")[1:])
            scaffold = key.split("[")[1].split(":")[0]
            start = key.split("[")[1].split(":")[1].split("-")[0]
            end = key.split("[")[1].split(":")[1].split("-")[1].split(" ")[0]
            strand = key.split("[")[1].split(" ")[1]
            kdomains = key.split("K=")[1]
            length = len(fa.full_dict[key])
            seq = fa.full_dict[key]
            out.write("%s\t%s\t%s\t%s\t%s\t'%s\t%s\t%s\t%s\n" %(id, description, scaffold, start, end,
                                                               strand, kdomains, length, seq))

def phylogram_fasta_mapping(fasta,out_label):
    fa = Fasta(fasta)
    count_dict = {"Known protein":0, "Ixodes hexagonus":0, "Ixodes pacificus":0, "Ixodes scapularis":0, "Ixodes persulcatus":0, "Ixodes ricinus":0,"Amblyomma maculatum":0}
    with open("%s_mapping.txt" %out_label, "wt") as out:
        with open("%s_renamed.fa" % out_label, "wt") as out2:
            out.write("Protein\tSpecie\n")
            for key in fa.full_dict:
                if "_KP" in key:
                    out.write("%s\tKnown protein\n" % key.replace(" ","_"))
                    out2.write(">%s\n%s\n" % (key.replace(" ","_"), fa.full_dict[key]))
                    count_dict["Known protein"] += 1
                elif key.startswith("Ihex"):
                    if ";" in key:
                        name = key.split(";")[1]
                    else:
                        name=key
                    if "Name=" in name:
                        name = name.replace("Name=","")
                    name = name.replace(" ","_").replace("Ihex","")
                    out.write("Ihex_%s\tIxodes hexagonus\n" %name)
                    out2.write(">Ihex_%s\n%s\n" % (name, fa.full_dict[key]))
                    count_dict["Ixodes hexagonus"] += 1
                elif key.startswith("Ipac"):
                    if ";" in key:
                        name = key.split(";")[1]
                    else:
                        name=key
                    if "Name=" in name:
                        name = name.replace("Name=", "")
                    name = name.replace(" ","_").replace("Ipac","")
                    out.write("Ipac_%s\tIxodes pacificus\n" %name)
                    out2.write(">Ipac_%s\n%s\n" % (name, fa.full_dict[key]))
                    count_dict["Ixodes pacificus"] += 1
                elif key.startswith("ISCP"):
                    out.write("Isca%s\tIxodes scapularis\n" %key.replace(" ","_").replace("ISCP",""))
                    out2.write(">Isca%s\n%s\n" % (key.replace(" ", "_").replace("ISCP", ""), fa.full_dict[key]))
                    count_dict["Ixodes scapularis"] += 1
                elif key.startswith("Iper"):
                    if ";" in key:
                        name = key.split(";")[1]
                    else:
                        name=key
                    if "Name=" in name:
                        name = name.replace("Name=","")
                    name = name.replace(" ","_").replace("Iper","")
                    out.write("Iper_%s\tIxodes persulcatus\n" %name)
                    out2.write(">Iper_%s\n%s\n" % (name, fa.full_dict[key]))
                    count_dict["Ixodes persulcatus"] += 1
                elif key.startswith("Iric"):
                    if not key.startswith("Iric_"):
                        name = key.replace(" ","_").replace("Iric","")
                    else:
                        name= key.replace(" ","_").replace("Iric_","")
                    out.write("Iric_%s\tIxodes ricinus\n" %name)
                    out2.write(">Iric_%s\n%s\n" % (name, fa.full_dict[key]))
                    count_dict["Ixodes ricinus"] += 1
                elif key.startswith("Amac"):
                    if not key.startswith("Amac_"):
                        name = key.replace(" ","_").replace("Amac","")
                    else:
                        name= key.replace(" ","_").replace("Amac_","")
                    out.write("Amac_%s\tAmblyomma maculatum\n" %name)
                    out2.write(">Amac_%s\n%s\n" % (name, fa.full_dict[key]))
                    count_dict["Amblyomma maculatum"] += 1
                elif key.startswith("Isca"):
                    out.write("Isca%s\tIxodes scapularis\n" % key.replace(" ", "_").replace("Isca", ""))
                    out2.write(">Isca%s\n%s\n" % (key.replace(" ", "_").replace("Isca", ""), fa.full_dict[key]))
                    count_dict["Ixodes scapularis"] += 1
                else:
                    out.write("Iric_%s\tIxodes ricinus\n" %key.replace(" ","_"))
                    out2.write(">Iric_%s\n%s\n" % (key.replace(" ", "_"), fa.full_dict[key]))
                    count_dict["Ixodes ricinus"] += 1
    with open("%s_count.txt" %out_label, "wt") as out:
        out.write("Species\tNumber of proteins\n")
        for key in count_dict:
            out.write("%s\t%i\n" %(key, count_dict[key]))

def getrow(line):
    row = line.replace("\n","").split("\t")
    return row

def get_miRNAs_from_gff(gff, out_label):
    with open(gff) as fl:
        gff = list(fl)

    count_dict = dict()
    for line in gff:
        if "hairpin" in line:
            row = line.replace("\n","").split("\t")
            mirna = row[8].split("gene_name ")[1].replace(";","")
            if mirna not in count_dict:
                count_dict[mirna] = 1
            else:
                count_dict[mirna] += 1

    with open("%s_miRNAs.txt" %out_label, "wt") as out1:
        with open("%s_miRNAs_numAlignments.tsv" %out_label, "wt") as out2:
            out2.write("miRNA\tNum Alignments\n")
            for key in count_dict:
                out1.write("%s\n" %key)
                out2.write("%s\t%i\n" %(key, count_dict[key]))

    with open("%s_final.gff" % out_label, "wt") as out:
        for line in gff:
            if "hairpin" in line:
                row = line.replace("\n", "").split("\t")
                scaff = row[0]
                software = row[1]
                type = row[2]
                start = row[3]
                end = row[4]
                score = row[5]
                strand = row[6]
                frame = row[7]
                description = row[8]
                mirna = row[8].split("gene_name ")[1].replace(";", "")
                if count_dict[mirna] < 3:
                    name = "%s_%s_%s" %(mirna, scaff, start)
                    out.write(
                        "%s\tsRNAbench\tgene\t%s\t%s\t%s\t%s\t%s\tID=%s;type=miRNA;\n" % (
                        scaff, start, end, score, strand,
                        frame, name))
                    out.write(
                        "%s\tsRNAbench\ttranscript\t%s\t%s\t%s\t%s\t%s\tID=%s.1;Parent=%s;type=miRNA;\n" % (
                            scaff, start, end, score, strand,
                            frame, name, name))
                    out.write(
                        "%s\tsRNAbench\texon\t%s\t%s\t%s\t%s\t%s\tID=%s.1.exon1;Parent=_%s.1;type=miRNA;\n" % (
                            scaff, start, end, score, strand,
                            frame, name, name))
            elif "mature" in line:
                if count_dict[mirna] < 3:
                    mature_row = line.replace("\n", "").split("\t")
                    start = mature_row[3]
                    end = mature_row[4]
                    mature_name = mature_row[8].split("gene_name ")[1].replace(";", "")
                    out.write(
                        "%s\tsRNAbench\tmature\t%s\t%s\t%s\t%s\t%s\tID=%s_%s_%s;Parent=_%s.1.exon1;type=miRNA_mature;\n" % (
                            scaff, start, end, score, strand,
                            frame, mature_name, scaff, start, name))

def extract_rRNA(fasta, output):
    fa = Fasta(fasta)
    with open(output, "wt") as out:
        for key in fa.full_dict:
            if "rRNA" in key:
                out.write(">%s\t%s\n" %(key, fa.full_dict[key]))

def depurate_known_fasta(fasta, output):
    fa = Fasta(fasta)
    with open(output, "wt") as out:
        for key in fa.full_dict:
            name = "_".join(key.split(",")[0].split(" ")[1::])
            out.write(">%s_KP\n%s\n" %(name, fa.full_dict[key]))

def signalp6_fa_depurate(signalp_fa, fasta, output):
    fa = Fasta(fasta)
    sigfa = Fasta(signalp_fa)

    with open(output, "wt") as out:
        for key in fa.full_dict:
            id = key.split(" ")[0]
            if id in sigfa.dict:
                out.write(">%s\n%s\n" %(key, sigfa.dict[id]))
            else:
                out.write(">%s\n%s\n" % (key, fa.dict[id]))

def kunitz_fasta_names(fasta, gff3, output):
    fa = Fasta(fasta)

    name_dict = dict()

    with open(gff3) as fl:
        gff = list(fl)

    k = 0
    for line in gff:
        k += 1
        if not line.startswith("#") and not line=="\n" and not line=="":
            row = line.split("\t")
            if row[2] == "mRNA":
                id = line.split("ID=")[1].split(";")[0]
                name = line.split("Name=")[1].split(";")[0]
                if " - Part " in name:
                    name = name.split(" - Part")[0]
                if "Partial" in name:
                    name = name.replace(" - Partial", "_partial")
                if name in name_dict.values():
                    for k in range(1,1000):
                        if "%s %i" %(name,k) not in name_dict.values():
                            name_dict[id] = "%s %i" %(name,k)
                            break
                else:
                    name_dict[id] = name
    with open(output,"wt") as out:
        for key in fa.full_dict:
            id = key.split(" ")[0]
            if id in name_dict:
                out.write(">%s %s\n%s\n" %("_".join(name_dict[id].split(" ")), key, fa.full_dict[key]))
            else:
                out.write(">Unnamed %s\n%s\n" % (key, fa.full_dict[key]))

def kunitz_by_domain_number(fasta, out_label):
    fa = Fasta(fasta)
    kn_dict = dict()
    for key in fa.full_dict:
        kn = "K%s" %key.split("K=")[1]
        if kn not in kn_dict:
            kn_dict[kn] = [key]
        else:
            kn_dict[kn].append(key)
    print(kn_dict)
    for key in kn_dict:
        with open("%s_%s.fa" %(out_label,key),"wt") as out:
            for id in kn_dict[key]:
                out.write(">%s\n%s\n" %(id, fa.full_dict[id]))

def kunitz_and_cystatins(fasta, family_file, domain_file, outlabel):
    fa = Fasta(fasta)
    family_kunitz = []
    family_cystatin = []
    domain_kunitz = []
    domain_cystatin = []

    with open(family_file) as ffile:
        family_lines = list(ffile)
    with open(domain_file) as dfile:
        domain_lines = list(dfile)

    for line in family_lines:
        if "Kunitz" in line or "kunitz" in line:
            family_kunitz.append(line.split("\t")[0])
        if "Cystatin" in line or "cystatin" in line:
            family_cystatin.append(line.split("\t")[0])

    for line in domain_lines:
        if "Kunitz" in line or "kunitz" in line:
            domain_kunitz.append(line.split("\t")[0])
        if "Cystatin" in line or "cystatin" in line:
            domain_cystatin.append(line.split("\t")[0])
    kunitz = deduplicate_list(family_kunitz + domain_kunitz)
    cystatin = deduplicate_list(family_cystatin + domain_cystatin)

    kunitz_outname = "%s_Kunitz.pep" %outlabel
    cystatin_outname = "%s_Cystatin.pep" %outlabel
    with open(kunitz_outname, "wt") as out_kunitz:
        with open(cystatin_outname, "wt") as out_cystatin:
            for id in fa.dict:
                if id in kunitz:
                    out_kunitz.write(">%s\n%s\n" %(id, fa.dict[id]))
                if id in cystatin:
                    out_cystatin.write(">%s\n%s\n" %(id, fa.dict[id]))

def defensins(fasta, family_file, domain_file, outlabel):
    fa = Fasta(fasta)
    family_defensin = []
    domain_defensin = []

    with open(family_file) as ffile:
        family_lines = list(ffile)
    with open(domain_file) as dfile:
        domain_lines = list(dfile)

    for line in family_lines:
        if "Defensin" in line or "denfensin" in line:
            family_defensin.append(line.split("\t")[0])

    for line in domain_lines:
        if "Defensin" in line or "defensin" in line:
            domain_defensin.append(line.split("\t")[0])
    defensin = deduplicate_list(family_defensin + domain_defensin)

    defensin_outname = "%s_Defensin.pep" %outlabel
    with open(defensin_outname, "wt") as out_defensin:
        for id in fa.dict:
            if id in defensin:
                out_defensin.write(">%s\n%s\n" %(id, fa.dict[id]))

def kunitz_annotate_ndomains(fasta, domain_file, output):

    ndomain_dict = dict()
    tmp_dict = dict()

    with open(domain_file) as dfile:
        domain_lines = list(dfile)

    for line in domain_lines:
        if "Kunitz" in line or "kunitz" in line:
            pep = line.split("\t")[0]
            program = line.split("\t")[3]
            if program != "PRINTS":
                if pep not in tmp_dict:
                    tmp_dict[pep] = [program]
                else:
                    tmp_dict[pep].append(program)

    for pep in tmp_dict:
        ndomain_dict[pep] = tmp_dict[pep].count(max(tmp_dict[pep], key=tmp_dict[pep].count))

    fa = Fasta(fasta)

    with open(output,"wt") as out:
        for key in fa.full_dict:
            name = key.split(" ")[0]
            des = "%s K=%i" %(key, ndomain_dict[name])
            out.write(">%s\n%s\n" %(des, fa.full_dict[key]))

def kunitz_ion_channel(fasta, blast, blastmap, output):
    fa = Fasta(fasta)
    blast = Outfmt6(blast)
    potassium_list = []
    potassium_peptides = []

    with open(blastmap) as bm:
        bm_lines = list(bm)
    for k in range(1, len(bm_lines)):
        if "potassium channel" in bm_lines[k] or "Potassium channel" in bm_lines[k]:
            id = bm_lines[k].split("\t")[0]
            potassium_list.append(id)
    for id in blast.pname_dict:
        if blast.pname_dict[id] in potassium_list:
            potassium_peptides.append(id)

    fa.filter(potassium_peptides, output)

# def BlastMapAnnotation(blast, blastmap, output)
#     blast = Outfmt6(blast)
#
#     with open(blastmap) as bm:
#         bm_lines = list(bm)
#
#     for id in blast.pname_dict:

def getncRNA_gffcompare(input_label, output_label, gff, fasta):
    coding_transcripts = []
    noncoding_transcripts = []
    if exists("%s=.tsv" %input_label):
        file = GffCompareTrackingGroup("%s=.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%sc.tsv" %input_label):
        file = GffCompareTrackingGroup("%sc.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%sk.tsv" %input_label):
        file = GffCompareTrackingGroup("%sk.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%sm.tsv" %input_label):
        file = GffCompareTrackingGroup("%sm.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%sn.tsv" %input_label):
        file = GffCompareTrackingGroup("%sn.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%sj.tsv" %input_label):
        file = GffCompareTrackingGroup("%sj.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%se.tsv" %input_label):
        file = GffCompareTrackingGroup("%se.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%so.tsv" %input_label):
        file = GffCompareTrackingGroup("%so.tsv" %input_label)
        coding_transcripts = coding_transcripts + file.transcripts
    if exists("%ss.tsv" %input_label):
        file = GffCompareTrackingGroup("%ss.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts
    if exists("%sx.tsv" %input_label):
        file = GffCompareTrackingGroup("%sx.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts
    if exists("%si.tsv" %input_label):
        file = GffCompareTrackingGroup("%si.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts
    if exists("%sy.tsv" %input_label):
        file = GffCompareTrackingGroup("%sy.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts
    if exists("%sp.tsv" %input_label):
        file = GffCompareTrackingGroup("%sp.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts
    if exists("%sr.tsv" %input_label):
        file = GffCompareTrackingGroup("%sr.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts
    if exists("%su.tsv" %input_label):
        file = GffCompareTrackingGroup("%su.tsv" %input_label)
        noncoding_transcripts = noncoding_transcripts + file.transcripts

    coding_transcripts = deduplicate_list(coding_transcripts)
    noncoding_transcripts = deduplicate_list(noncoding_transcripts)

    print("Conding transcripts = %i" %len(coding_transcripts))
    print("Non conding transcripts = %i" % len(noncoding_transcripts))

    noncoding_transcripts = [x for x in noncoding_transcripts if x not in coding_transcripts]

    print("Exclusive Non conding transcripts = %i" % len(noncoding_transcripts))

    #Writing noncoding transcripts
    with open("%s.txt" %output_label, "wt") as out:
        out.write("\n".join(noncoding_transcripts))

    noncoding_transcripts_tofilter = ["_".join(x.split("_")[0:-1]) for x in noncoding_transcripts]

    noncoding_transcripts = [x.replace("_mid", ".aln") for x in noncoding_transcripts]
    ## Filtering gff.
    with open(gff) as input_gff:
        gff_lines = list(input_gff)
    with open("%s.gff" %output_label, "wt") as out:
        out.write(gff_lines[0])
        for k in range(1,len(gff_lines)):
            if "ID="  in gff_lines[k]:
                id = gff_lines[k].split("ID=")[1].split(";")[0]
            if ".exon" in id:
                id = ".".join(id.split(".")[0:-1])
            if id in noncoding_transcripts:
                out.write(gff_lines[k])

    ## Filtering Fasta:
    fa = Fasta(fasta)
    fa.filter(noncoding_transcripts, "%s.fa" %output_label)

def stringDB_enrichment_figure(bp_file,mf_file,cc_file, cl_file, out):
    description_list =[]
    fdr_list = []
    strength_list=[]
    if bp_file:
        bp_fl = StringDBGO(bp_file)
        for key in bp_fl.fdr_dict:
            description_list.append(key)
            fdr_list.append(bp_fl.fdr_dict[key])
            strength_list.append(bp_fl.strength_dict[key])
    if mf_file:
        mf_fl = StringDBGO(mf_file)
        for key in mf_fl.fdr_dict:
            description_list.append(key)
            fdr_list.append(mf_fl.fdr_dict[key])
            strength_list.append(mf_fl.strength_dict[key])
    if cc_file:
        cc_fl = StringDBGO(cc_file)
        for key in cc_fl.fdr_dict:
            description_list.append(key)
            fdr_list.append(cc_fl.fdr_dict[key])
            strength_list.append(cc_fl.strength_dict[key])
    if cl_file:
        cl_fl = StringDBGO(cl_file)
        for key in cl_fl.fdr_dict:
            description_list.append(key)
            fdr_list.append(cl_fl.fdr_dict[key])
            strength_list.append(cl_fl.strength_dict[key])
        

    df = pd.DataFrame({'-log10(FDR)': fdr_list, 'Strength': strength_list}, index=description_list)
    df = df.reindex(index=df.index[::-1])

    sns.set_theme()
    sns.set_style("dark")
    fig = plt.figure()  # Create matplotlib figure

    ax = fig.add_subplot(111)  # Create matplotlib axes
    ax2 = ax.twiny()  # Create another axes that shares the same x-axis as ax.

    width = 0.4

    df["-log10(FDR)"].plot.barh(color="royalblue", ax=ax, width=width, position=1, fontsize="18")
    df["Strength"].plot.barh(color="darkorange", ax=ax2, width=width, position=0, fontsize="18")

    ax.set_xlabel('-log10(FDR)', fontsize = "18")
    ax2.set_xlabel('Strength', fontsize = "18")

    #legend
    custom_lines = [Line2D([0], [0], color="royalblue", lw=4),
                    Line2D([0], [0], color="darkorange", lw=4)]
    ax2.legend(custom_lines, ['-log10(FDR)', 'Strength'], loc="upper right", fontsize="16")
    #ax = df.plot.barh()
    #fig = ax.get_figure()
    fig.set_figheight(10)
    fig.set_figwidth(15)
    plt.tight_layout()
    fig.savefig(out)
    plt.close()

def intersect_masigpro_clusters(file1, file2, out_label):
    with open(file1) as fl1:
        fl1 = list(fl1)
    fl1 = [x.replace("\n","") for x in fl1]
    with open(file2) as fl2:
        fl2 = list(fl2)
    fl2 = [x.replace("\n", "") for x in fl2]
    intersected_list = intersection(fl1, fl2)
    not_intersected_list_mg = intersection(fl1, fl2, True)
    not_intersected_list_sg = intersection(fl2, fl1, True)
    print("Intersected DEGs: %i" %len(intersected_list))
    print("Exclusive DEGs for MG: %i" %len(not_intersected_list_mg))
    print("Exclusive DEGs for SG: %i" % len(not_intersected_list_sg))

    with open("%s_intersection.txt" %out_label, "wt") as out:
        out.write("\n".join(intersected_list))
    with open("%s_exclusive_MG.txt" %out_label, "wt") as out:
        out.write("\n".join(not_intersected_list_mg))
    with open("%s_exclusive_SG.txt" %out_label, "wt") as out:
        out.write("\n".join(not_intersected_list_sg))

    plt.figure(figsize=(5, 5))
    venn2(subsets=(len(not_intersected_list_mg), len(not_intersected_list_sg), len(intersected_list)),
          set_labels=("MG DEGs", "SG DEGs"))
    plt.savefig("%s_venn.png" % out_label)
    plt.close()

def distribution_by_clusters(file, cluster_label, output):
    with open(file) as fl:
        fl = list(fl)
    ids = [x.replace("\n","") for x in fl]
    cl_dict = dict()
    for k in range(1,10):
        with open("%s%i.txt" %(cluster_label, k)) as cl_file:
            cl_file = list(cl_file)
        cl_list = [x.replace("\n","") for x in cl_file]
        cl_dict["Cluster %i" %k] = 0
        for id in cl_list:
            if id in ids:
                cl_dict["Cluster %i" % k] += 1
        cl_dict["Cluster %i" % k] += cl_dict["Cluster %i" % k] / len(cl_list) *100

    for key in cl_dict:
        print("%s: %i" %(key, cl_dict[key]))

    barplot_from_dict(cl_dict, output, "Percentage of DEGs common in MG and SG")

def intersectfasta(fasta, excel, output, reverse=False):
    ft = Fasta(fasta)
    ex = Excel(excel)
    with open(output, "wt") as out:
        for key in ft.dict:
            if not reverse:
                if key in ex.IDs:
                    out.write(">%s\n%s\n" % (key, ft.dict[key]))
            elif reverse:
                if key not in ex.IDs:
                    out.write(">%s\n%s\n" % (key, ft.dict[key]))

def intersectexcel(excel, fasta, output, reverse=False):
    ft = Fasta(fasta)
    ex = CompleteExcel(excel)
    with open(output, "wt") as out:
        out.write(ex.lines[0])
        for n in range(1,len(ex.lines)):
            line = CompleteExcel.Line(ex.lines[n])
            if not reverse:
                if line.cdsID in ft.dict:
                    out.write(ex.lines[n])
            elif reverse:
                if line.cdsID not in ft.dict:
                    out.write(ex.lines[n])

def intersect_background(IDs,background,output):
    with open(background) as bck:
        bck = list(bck)
    with open(output,"wt") as out:
        for line in bck:
            ID = line.split("\t")[0].split(".")[0]
            if ID in IDs:
                out.write(line)

def gethitID(row):
    db = row[7]
    if db == "Uniref90":
        hit = row[2].split("_")[-1]
    elif db == "TickSialoFam":
        hit = row[2]
    elif db == "ArachnidaDB" or "SwissProt":
        hit = row[2].split("|")[-1]
    return hit

def get_expression_dicts(lines):
    evalues_secreted = []
    evalues_nonsecreted = []
    for n in range(1,len(lines)):
        line = CompleteExcel.Line(lines[n])
        if line.classification == "Putatively Secreted":
            evalues_secreted.append(float(line.total_FPKM))
        else:
            evalues_nonsecreted.append(float(line.total_FPKM))
    return evalues_secreted, evalues_nonsecreted

def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.

    Args:
        list_of_lengths (list): List of numbers.

    Returns:
        float: N50 value.

    """
    tmp = []
    for tmp_number in set(list_of_lengths):
        tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()

    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]

    return median

def writenames(names, output):
    with open(output, "wt") as out:
        out.write("\n".join(names))

def deduplicate_list(mylist):
    mylist = list(dict.fromkeys(mylist))
    return mylist

def getlength(dictionary):
    length_dict = dict()
    for key in dictionary:
        length_dict[key] = len(dictionary[key])
    return length_dict

def getUnknown(excel, uniprot, output):
    ex = Excel(excel)
    uni = Uniprot(uniprot)
    uni_IDs = uni.getAnnotated()
    with open(output, "wt") as out:
        for n in range(1, len(ex.lines)):
            line = Excel.Line(ex.lines[n])
            if line.hitID not in uni_IDs:
                out.write(">%s\n%s\n" % (line.cdsID, line.seq))

def startbymethionine(fasta, output):
    ft = Fasta(fasta)
    count = 0
    total_count = 0
    with open(output, "wt") as out:
        for key in ft.dict:
            total_count += 1
            if ft.dict[key].startswith("M"):
                count += 1
                out.write(">%s\n%s\n" % (key, ft.dict[key]))
    print("%s of a total of %s CDS start with methionine." % (count, total_count))

def mergegoterms(goterm1, goterm2):
    merged = goterm1.split("; ") + goterm2.split("; ")
    dd = deduplicate_list(merged)
    return dd

def createCompleteExcel(excel, filtered_fasta, cds_fasta, cds_pep, matrix_FPKM, uniprot, blast2go, interpro, signalp, phobius, tmhmm_mature, tmhmm, output, GOobo_file, list_of_DEfiles):
    #Loading ildes
    ex = Excel(excel)
    ft = Fasta(filtered_fasta)
    pep = Fasta(cds_pep)
    cds = Fasta(cds_fasta)
    fpkm = Ematrix(matrix_FPKM)
    obo = GOobo(GOobo_file)
    uni = Uniprot(uniprot, obo.dict)
    b2go = Blast2GO(blast2go, obo.dict)
    ip = Interpro(interpro, obo.dict)
    tm_mature = TMHMM(tmhmm_mature)
    sp = SignalP4(signalp, tm_mature)
    tm = TMHMM(tmhmm)
    ph = Phobius(phobius)
    DE_dict, FC_dict, pvalue_dict = DEsummarize(list_of_DEfiles)
    uni_count = 0
    b2go_count = 0

    # Header
    with open(output, "wt") as out:
        out.write("ID\tTranscript Sequence\tCDS Sequence\tPeptide Sequence\tFamily\tBestHitID\tBestHitIdentity\tBestHitCoverage\tBestHitBitScore\tBestHitDescription\tBestHitDatabase\tUniprotGO\tUniprotGO names\tKeywords\tKeywords names\tInterproID\tInterpro names\tInterproGO\tInterproGO names\tMergedGO\tMergedGO names\tSignalP\tTMHMM in mature (only for SP)\tTMHMM\tClass\tTissue Specifity\tSG_FC\tMG_FC\tTotal FPKMs\tDE analysis\tDE FC\tPosterior probability of being DE (PPDE)\t%s\n" % ("\t".join(fpkm.order)))
        for key in cds.dict:
            cdsID = key
            transcriptID = key.split(".")[0]
            cds_seq = cds.dict[key]
            pep_seq = pep.dict[key]
            transcript_seq = ft.dict[transcriptID]

            # Excel
            if key in ex.dict:
                ex_line = Excel.Line(ex.dict[key])
                hitID = gethitID(ex_line.row)
                identity = ex_line.row[3]
                coverage = ex_line.row[4]
                bitscore = ex_line.row[5]
                description = ex_line.row[6]
                db = ex_line.row[7]
            else:
                hitID = "None"
                identity = "None"
                coverage = "None"
                bitscore = "None"
                description = "None"
                db = "None"

            if hitID == "V5H0G4":
                print(cdsID)
                print(uni.goterm_dict[hitID])

            # Uniprot and Blast2GO
            if hitID in uni.goterm_dict and uni.goterm_dict[hitID] != "None" and uni.goterm_dict[hitID]!= "":
                GOterm = uni.goterm_dict[hitID]
                GOterm_des = uni.goterm_des_dict[hitID]
                uni_count += 1
            elif cdsID in b2go.dict and b2go.dict[cdsID] != "no GO terms" and b2go.dict[cdsID] != "":
                GOterm = b2go.dict[cdsID]
                GOterm_des = b2go.des_dict[cdsID]
                b2go_count += 1
            else:
                GOterm = "None"
                GOterm_des = "None"
            if hitID in uni.keyword_dict and uni.keyword_dict[hitID] != "None" and uni.keyword_dict[hitID]!= "":
                KW = uni.keyword_dict[hitID]
                KW_des = uni.keyword_des_dict[hitID]
            else:
                KW = "None"
                KW_des = "None"

            # Interpro
            if cdsID in ip.ipro_dict:
                ipro_id = ip.ipro_dict[cdsID]
                ipro_des = ip.ipro_des_dict[cdsID]
            else:
                ipro_id = "None"
                ipro_des = "None"
            if cdsID in ip.goterm_dict:
                ipro_goterm = ip.goterm_dict[cdsID]
                ipro_goterm_des = ip.goterm_des_dict[cdsID]
            else:
                ipro_goterm = "None"
                ipro_goterm_des = "None"

            # Merge GO terms
            if ipro_goterm != "None" and GOterm != "None":
                merged = "; ".join(mergegoterms(GOterm, ipro_goterm))
                merged_des = "; ".join(mergegoterms(GOterm_des, ipro_goterm_des))
            elif ipro_goterm == "None" and GOterm != "None":
                merged = GOterm
                merged_des = GOterm_des
            elif ipro_goterm != "None" and GOterm == "None":
                merged = ipro_goterm
                merged_des = ipro_goterm_des
            else:
                merged = "None"
                merged_des = "None"

            #Family
            if db == "Uniref90" or db == "SwissProt" or db == "ArachnidaDB":
                if hitID in uni.family_dict and uni.family_dict[hitID] != "None":
                    family = uni.family_dict[hitID]
                else:
                    family = " ".join(description.split("=")[0].split(" ")[1:-1])
            elif db == "TickSialoFam":
                family = "%s family" %description.split("|")[0].split(", ")[1]
            elif db == "None" and ipro_des != "None":
                family = ipro_des.split("; ")[0]
            elif db == "None" and merged_des != "None":
                family = merged_des.split("; ")[0]
            else:
                family = "None"



            # SignalP and Phobius
            if cdsID in sp.dict:
                sigp = sp.dict[cdsID]
                if cdsID in tm_mature.dict:
                    sigp_mature_tmhmm = tm_mature.dict[cdsID]
                else:
                    sigp_mature_tmhmm = "None"
            else:
                sigp = "None"
                sigp_mature_tmhmm = "None"

            if cdsID in ph.dict:
                if ph.dict[cdsID] == "Yes" and sigp == "None":
                    sigp = ph.dict[cdsID]


            # TMHMM
            if cdsID in tm.dict:
                tm_domain = tm.dict[cdsID]
            else:
                tm_domain = "None"

            # Class
            if sigp != "None" and sigp_mature_tmhmm == "None" and (merged != "None" or ipro_id != "None"):
                classification = "Putatively Secreted, Annotated"
            elif sigp != "None" and sigp_mature_tmhmm == "None" and merged == "None" and ipro_id == "None":
                classification = "Putatively Secreted, not Annotated"
            elif sigp != "None" and sigp_mature_tmhmm != "None" and merged == "None" and ipro_id == "None":
                classification = "Putatively Non-Secreted, not Annotated"
            elif sigp != "None" and sigp_mature_tmhmm != "None" and (merged != "None" or ipro_id != "None"):
                classification = "Putatively Non-Secreted, Annotated"
            elif sigp == "None" and (merged != "None" or ipro_id != "None"):
                classification = "Putatively Non-Secreted, Annotated"
            elif sigp == "None" and merged == "None" and ipro_id == "None":
                classification = "Putatively Non-Secreted, not Annotated"

            # Modificating classification according to the blast hit
            if "Secreted" in description or "secreted" in description or "|S|" in description:
                if merged != "None" or ipro_id != "None":
                    classification = "Putatively Secreted, Annotated"
                else:
                    classification = "Putatively Secreted, not Annotated"

            # Expression
            evalues = fpkm.dict[transcriptID]
            # Expression class
            evalues_float = []
            for item in evalues:
                evalues_float.append(float(item))
            total_FPKM = np.sum(evalues_float)
            SG_evalues = evalues_float[0:55]
            MG_evalues = evalues_float[55:88]
            SG_mean = np.mean(SG_evalues)
            MG_mean = np.mean(MG_evalues)
            SG_ratio = SG_mean / (MG_mean + SG_mean)
            MG_ratio = 1 - SG_ratio
            SG_FC = np.log2(SG_mean/MG_mean)
            MG_FC = np.log2(MG_mean/SG_mean)
            if math.isinf(SG_FC) and SG_FC < 0:
                SG_FC = -9.99999999999999*(10**307)
            elif math.isinf(SG_FC) and SG_FC > 0:
                SG_FC = 9.99999999999999*(10**307)
            if math.isinf(MG_FC) and MG_FC < 0:
                MG_FC = -9.99999999999999*(10**307)
            elif math.isinf(MG_FC) and MG_FC > 0:
                MG_FC = 9.99999999999999*(10**307)
            if SG_FC >= 2:
                eclass = "SG specific"
            elif MG_FC >= 2:
                eclass = "MG specific"
            else:
                eclass = "Neutral"

            #Differential Expression
            if transcriptID in DE_dict:
                DE_str = "; ".join(DE_dict[transcriptID])
                FC = "; ".join(FC_dict[transcriptID])
                pvalue = "; ".join(pvalue_dict[transcriptID])
            else:
                DE_str = "Stable"
                FC = "None"
                pvalue = "None"



            # Writing
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cdsID, transcript_seq, cds_seq, pep_seq, family, hitID, identity, coverage, bitscore, description, db, GOterm, GOterm_des, KW, KW_des, ipro_id, ipro_des, ipro_goterm, ipro_goterm_des, merged, merged_des, sigp, sigp_mature_tmhmm, tm_domain, classification, eclass, float(SG_FC), float(MG_FC), total_FPKM, DE_str, FC, pvalue, "\t".join(evalues)))
    print ("Uniprot GO terms: %s\nB2GO GO terms: %s\n" %(uni_count, b2go_count))

def DEsummarize(list_of_files):
    # The file must have two columns, one for the name of the file, and one for the coding it will have in the string.
    DE_dict = {}
    FC_dict = {}
    pvalue_dict = {}

    with open(list_of_files) as filelist:
        filelist = list(filelist)
        for line in filelist:
            file = line.split("\t")[0]
            name = line.split("\t")[1].replace("\n","")
            if "MG" in file.split("/")[-1]:
                tissue = "MG"
            elif "SG" in file.split("/")[-1]:
                tissue = "SG"
            de = DE(file)
            for key in de.IDsUPfc:
                if key not in DE_dict:
                    DE_dict[key] = ["%s|%sUP" %(tissue,name)]
                    FC_dict[key] = [str(de.IDsUPfc[key])]
                    pvalue_dict[key] = [str(de.IDsUPpvalue[key])]
                else:
                    DE_dict[key] += ["%s|%sUP" %(tissue,name)]
                    FC_dict[key] += [str(de.IDsUPfc[key])]
                    pvalue_dict[key] += [str(de.IDsUPpvalue[key])]
            for key in de.IDsDOWNfc:
                if key not in DE_dict:
                    DE_dict[key] = ["%s|%sDOWN" %(tissue,name)]
                    FC_dict[key] = [str(de.IDsDOWNfc[key])]
                    pvalue_dict[key] = [str(de.IDsDOWNpvalue[key])]
                else:
                    DE_dict[key] += ["%s|%sDOWN" %(tissue,name)]
                    FC_dict[key] += [str(de.IDsDOWNfc[key])]
                    pvalue_dict[key] += [str(de.IDsDOWNpvalue[key])]
    return DE_dict, FC_dict, pvalue_dict

def DEgetlist(list_of_files, excel_transcriptIDs=False):
    all_dict = {}
    with open(list_of_files) as filelist:
        filelist = list(filelist)
        for line in filelist:
            file = line.split("\t")[0]
            name = line.split("\t")[1].replace("\n","").replace("|","-")
            de = DE(file, excel_transcriptIDs)
            UP_name = "%sUP.txt" %name
            DOWN_name = "%sDOWN.txt" %name
            if UP_name not in all_dict:
                all_dict[UP_name] = de.IDsUP
            else:
                all_dict[UP_name] += de.IDsUP
            if DOWN_name not in all_dict:
                all_dict[DOWN_name] = de.IDsDOWN
            else:
                all_dict[DOWN_name] += de.IDsDOWN
    for key in all_dict:
        with open(key,"wt") as out:
            out.write("\n".join(all_dict[key]))
    files = list(all_dict.keys())
    return files

def numberofB2GO(excel,b2go,obo):
    ex = CompleteExcel(excel)
    obo = GOobo(obo)
    b2go = Blast2GO(b2go, obo.dict)
    count = 0
    for n in range(1,len(ex.lines)):
        line = CompleteExcel.Line(ex.lines[n])
        if line.cdsID in b2go.dict and b2go.dict[line.cdsID] != "None":
            count +=1
    print("Number of B2GO hits: %s" %count)

def enrichment(files,reference_label, pwd, comparison=False):
    os.system("mkdir %s/GO %s/KW %s/IP %s/family" %(pwd,pwd,pwd,pwd))
    for file in files:
        outname = "%s_goTerm.tsv" %file.split(".")[0]
        if comparison != False:
            output_label = "_".join(file.split("_")[0:-1])
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_GO_all.txt" %(file, output_label))
        os.system("mv %s GO" %outname)
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_KW_all.txt" %(file, output_label))
        os.system("mv %s KW" %outname)
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_IP_all.txt" %(file, output_label))
        os.system("mv %s IP" %outname)
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_family_all.txt" %(file, output_label))
        os.system("mv %s family" %outname)

def enrichment_command(file,reference, pwd, annotation):
    os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s" %(file, reference))
    outname = "%s_goTerm.tsv" %file.split(".")[0]
    folder = "%s/%s" %(pwd,annotation)

    os.system("mv %s %s/" %(outname, folder))

def enrichment_output_analysis(list_of_files, output_label):
    GO_enriched = {}
    GO_depleted = {}
    with open(list_of_files) as list_ft:
        list_ft = list(list_ft)
    for file in list_ft:
        file = file.replace("\n","")
        name = file.split("/")[-1].replace(".tsv","").replace("_goTerm","")
        enrichment_out = EnrichmentOut(file)
        GO_dict = enrichment_out.getGO()
        for GO in GO_dict:
            if GO_dict[GO][1] == "enriched":
                if GO in GO_enriched:
                    GO_enriched[GO][1] = GO_enriched[GO][1] + "; " + name
                    GO_enriched[GO][2] = GO_enriched[GO][2] + 1
                else:
                    GO_enriched[GO]=[GO_dict[GO][0],name,1]
            elif GO_dict[GO][1] == "depleted":
                if GO in GO_depleted:
                    GO_depleted[GO][1] = GO_depleted[GO][1] + "; " + name
                    GO_depleted[GO][2] = GO_depleted[GO][2] + 1
                else:
                    GO_depleted[GO]=[GO_dict[GO][0],name,1]
    with open("%s_enriched.tsv" %output_label, "wt") as out:
        out.write("ID\tdescription\tgroups\tcount\tMG_count\tSG_count\n")
        for GO in GO_enriched:
            mg_count=0
            sg_count=0
            for condition in GO_enriched[GO][1].split(";"):
                if "MG" in condition:
                    mg_count+=1
                elif "SG" in condition:
                    sg_count+=1

            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(GO, GO_enriched[GO][0], GO_enriched[GO][1], GO_enriched[GO][2],mg_count,sg_count))
    with open("%s_depleted.tsv" %output_label, "wt") as out:
        out.write("ID\tdescription\tgroups\tcount\tMG_count\tSG_count\n")
        for GO in GO_depleted:
            mg_count=0
            sg_count=0
            for condition in GO_depleted[GO][1].split(";"):
                if "MG" in condition:
                    mg_count+=1
                elif "SG" in condition:
                    sg_count+=1
            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(GO, GO_depleted[GO][0], GO_depleted[GO][1], GO_depleted[GO][2], mg_count, sg_count))

def gene_overlapping(list_of_files, output):
    overlap_dict = {}
    overlap_dict_group = {}
    with open(list_of_files) as flist:
        flist = list(flist)
    for file in flist:
        file = file.replace("\n","")
        name = file.split(".")[0]
        with open(file) as ft:
            ft = list(ft)
        for line in ft:
            ID = line.split("\t")[0].replace("\n","")
            if ID not in overlap_dict:
                overlap_dict_group[ID] = [name]
                overlap_dict[ID] = 1
            else:
                overlap_dict_group[ID].append(name)
                overlap_dict[ID] += 1
    overlap_dict = {k: v for k, v in sorted(overlap_dict.items(), key=lambda item: item[1], reverse=True)}
    with open(output,"wt") as out:
        out.write("ID\tOverlap count\tMG count\tSG count\tgroups\n")
        for key in overlap_dict:
            if overlap_dict[key] > 1:
                mg_count=0
                sg_count=0
                for item in overlap_dict_group[key]:
                    if "MG" in item:
                        mg_count += 1
                    if "SG" in item:
                        sg_count +=1
                out.write("%s\t%s\t%s\t%s\t%s\n" %(key, overlap_dict[key], mg_count, sg_count, overlap_dict_group[key]))

def spearman_corr_orderedbyweight(completeexcel, output):
    cexcel = CompleteExcel(completeexcel)
    MG_1_order, MG_2_order, SG_1_order, SG_2_order = cexcel.orderbyweight()
    MG_1_df = cexcel.df[MG_1_order]
    MG_2_df = cexcel.df[MG_2_order]
    SG_1_df = cexcel.df[SG_1_order]
    SG_2_df = cexcel.df[SG_2_order]
    MG_corr_dict = {}
    SG_corr_dict = {}

    for index, row in MG_1_df.iterrows():
        ID = index
        evalues_1 = list(row)
        evalues_2 = list(MG_2_df.loc[ID, :])
        sp = spearmanr(evalues_1, evalues_2)
        if not math.isnan(float(sp[0])):
            MG_corr_dict[ID] = sp
    for index, row in SG_1_df.iterrows():
        ID = index
        evalues_1 = list(row)
        evalues_2 = list(SG_2_df.loc[ID, :])
        sp = spearmanr(evalues_1, evalues_2)
        if not math.isnan(float(sp[0])):
            SG_corr_dict[ID] = sp

    #Ordering by correlation
    MG_corr_dict = {k: v for k, v in sorted(MG_corr_dict.items(), key=lambda item: item[1], reverse=True)}
    SG_corr_dict = {k: v for k, v in sorted(SG_corr_dict.items(), key=lambda item: item[1], reverse=True)}

    #Writting the output
    with open("%s_MG.tsv" %output, "wt") as out:
        out.write("ID\tSpearman Correlation\tP value\n")
        for key in MG_corr_dict:
            out.write("%s\t%f\t%f\n" %(key, MG_corr_dict[key][0], MG_corr_dict[key][1]))
    with open("%s_SG.tsv" %output, "wt") as out:
        out.write("ID\tSpearman Correlation\tP value\n")
        for key in SG_corr_dict:
            out.write("%s\t%f\t%f\n" %(key, SG_corr_dict[key][0], SG_corr_dict[key][1]))

def spearman_corr_byweight_distribution(file, output_label, plot_label, nbins=20):
    df = pd.read_csv(file, sep="\t")
    Plotting.plot_distribution(df,"Spearman Correlation","%s.pdf" %output_label, plot_label, nbins)

def spearman_corr_byweight_top(file, completeexcel, tissue, output_label, n=200):
    sp_dict = {}
    sp_top_dict = {}
    sp_bot_dict = {}
    with open(file) as sp:
        sp = list(sp)
    for k in range(1, len(sp)):
        ID = sp[k].split("\t")[0]
        sp_corr = sp[k].split("\t")[1]
        sp_dict[ID] = sp_corr

    cexcel = CompleteExcel(completeexcel)
    count = 0
    gcount = 0
    for key in sp_dict:
        gcount += 1
        if isDEtimecourse(cexcel.DE_dict[key],tissue):
            count += 1
            sp_top_dict[key] = [sp_dict[key], cexcel.DE_dict[key]]
        if gcount == n:
            nDE_top = count
        if count == n:
            break

    sp_dict_reverse = {k: v for k, v in sorted(sp_dict.items(), key=lambda item: item[1])}
    count = 0
    gcount = 0
    for key in sp_dict_reverse:
        gcount += 1
        if isDEtimecourse(cexcel.DE_dict[key], tissue):
            count += 1
            sp_bot_dict[key] = [sp_dict_reverse[key], cexcel.DE_dict[key]]
        if gcount == n:
            nDE_bot = count
        if count == n:
            break

    with open("%s_top%i.txt" %(output_label, n), "wt") as out:
        with open("%s_top%i_IDs.txt" % (output_label, n), "wt") as outID:
            out.write("ID\tSpearman Correlation\tDE analysis\n")
            for key in sp_top_dict:
                ID = key.split(".")[0]
                out.write("%s\t%s\t%s\n" %(key, sp_top_dict[key][0], "; ".join(sp_top_dict[key][1])))
                outID.write("%s\n" %ID)
    with open("%s_bot%i.txt" %(output_label, n), "wt") as out:
        with open("%s_top%i_IDs.txt" % (output_label, n), "wt") as outID:
            out.write("ID\tSpearman Correlation\tDE analysis\n")
            for key in sp_bot_dict:
                ID = key.split(".")[0]
                out.write("%s\t%s\t%s\n" %(key, sp_bot_dict[key][0], "; ".join(sp_bot_dict[key][1])))
                outID.write("%s\n" % ID)
    with open("%s_stats%i.txt" %(output_label,n), "wt") as out:
        out.write("Number of timecourse-DE genes in the top %i = %i\n" %(n,nDE_top))
        out.write("Number of timecourse-DE genes in the bot %i = %i\n" % (n, nDE_bot))

def isDEtimecourse(DE_list, tissue):
    isDE = False
    for item in DE_list:
        if item != "Stable":
            DEtissue = item.split("|")[0]
            exp1 = item.split("|")[1].split("_")[0]
            exp2 = item.split("|")[2].split("_")[0]
            if tissue == DEtissue and exp1 == exp2:
                isDE = True
    return isDE

def family_distribution_stats(cexcel, output):
    cexcel = CompleteExcel(cexcel)
    family_dict = {}
    family_perc_dict = {}
    family_n = 0
    for n in range(1, len(cexcel.lines)):
        line = CompleteExcel.Line(cexcel.lines[n])
        family = line.family.split(",")[0].replace("(","").replace(")","").replace("Putative","").replace("putative","").replace("Fragment","").replace("fragment","").replace("Conserved","").replace("conserved","").strip().capitalize()
        if family != "Uncharacterized protein" and family != "None":
            family_n += 1
            if family not in family_dict:
                family_dict[family] = [line.transcriptID]
            else:
                family_dict[family].append(line.transcriptID)

    for key in family_dict:
        family_perc_dict[key] = [len(family_dict[key])/family_n * 100, ",".join(family_dict[key])]
    print(family_n)
    df = pd.DataFrame.from_dict(family_perc_dict, orient="index", columns=["Family percentage","Transcript ID"])
    df = df.sort_values("Family percentage", ascending=False)
    df.to_csv(output, sep="\t")

def stability(ematrix, output, out_label, n_filter, m_filter, randomization=False):
    em = Ematrix(ematrix)
    index = []
    cv_dict = dict()
    cv_tuple_dict = dict()
    cv_list_sg = []
    expression_list_sg = []
    condition_list_sg = []
    exposure_list_sg = []
    transcript_list_sg = []
    cv_list_mg = []
    expression_list_mg = []
    condition_list_mg = []
    exposure_list_mg = []
    transcript_list_mg = []
    df_cv_ecdf_sg = pd.DataFrame()
    df_cv_ecdf_mg = pd.DataFrame()
    mg_cv_unstable_dict = dict()
    sg_cv_unstable_dict = dict()
    mg_rank_df = pd.DataFrame()
    mg_unfed_rank_df = pd.DataFrame()
    mg_1_rank_df = pd.DataFrame()
    mg_2_rank_df = pd.DataFrame()
    sg_rank_df = pd.DataFrame()
    sg_unfed_rank_df = pd.DataFrame()
    sg_1_rank_df = pd.DataFrame()
    sg_2_rank_df = pd.DataFrame()
    rank_df = pd.DataFrame()
    conditions = ["SG unfed", "SG first 12h", "SG first 24h", "SG first 48h", "SG first 72h", "SG first 96h", "SG second 12h", "SG second 24h",
                  "SG second 48h", "SG second 72h", "SG second 96h", "MG unfed", "MG first 12h", "MG first 24h", "MG first 48h", "MG first 72h", 
                  "MG first 96h", "MG second 12h", "MG second 24h", "MG second 48h", "MG second 72h", "MG second 96h"]
    SG_conditions = ["SG unfed", "SG first 12h", "SG first 24h", "SG first 48h", "SG first 72h", "SG first 96h", "SG second 12h", "SG second 24h",
                  "SG second 48h", "SG second 72h", "SG second 96h"]
    MG_conditions = ["MG unfed", "MG first 12h", "MG first 24h", "MG first 48h", "MG first 72h",
                  "MG first 96h", "MG second 12h", "MG second 24h", "MG second 48h", "MG second 72h", "MG second 96h"]
    for condition in conditions:
        cv_dict[condition] = []
        cv_tuple_dict[condition] = []
    for condition in SG_conditions:
        sg_cv_unstable_dict[condition] = []
    for condition in MG_conditions:
        mg_cv_unstable_dict[condition] = []
    for j in range(1,len(em.lines)):
        line = em.lines[j].replace("\n","").split("\t")
        ID = line[0]
        index.append(ID)
        for k in range(0,11):
            condition_evalues = line[k*5+3:k*5+6] #5+6 for 5 samples
            condition_evalues = [float(x) for x in condition_evalues]
            condition_list_sg.append(conditions[k])
            exposure_list_sg.append(conditions[k].split(" ")[1])
            transcript_list_sg.append(ID)
            if any([evalue>=5 for evalue in condition_evalues]):
                cv = variation(condition_evalues)
                cv_dict[conditions[k]].append(cv)
                cv_list_sg.append(cv)
                expression_list_sg.append(np.mean(condition_evalues))
                cv_tuple_dict[conditions[k]].append((ID,cv))
                if cv > n_filter:
                    sg_cv_unstable_dict[conditions[k]].append(ID)
            else:
                cv_dict[conditions[k]].append(NaN)
                cv_list_sg.append(NaN)
                expression_list_sg.append(NaN)
        for k in range(0,11):
            condition_evalues = line[k * 3 + 56:k * 3 + 59]
            condition_evalues = [float(x) for x in condition_evalues]
            condition_list_mg.append(conditions[k+11])
            exposure_list_mg.append(conditions[k+11].split(" ")[1])
            transcript_list_mg.append(ID)
            if any([evalue >= 5 for evalue in condition_evalues]):
                cv = variation(condition_evalues)
                cv_dict[conditions[k+11]].append(cv)
                cv_list_mg.append(cv)
                expression_list_mg.append(np.mean(condition_evalues))
                cv_tuple_dict[conditions[k+11]].append((ID, cv))
                if cv > n_filter:
                    mg_cv_unstable_dict[conditions[k+11]].append(ID)

            else:
                cv_dict[conditions[k + 11]].append(NaN)
                cv_list_mg.append(NaN)
                expression_list_mg.append(NaN)

    total_sg_genes = sg_cv_unstable_dict[conditions[0]]
    total_mg_genes = mg_cv_unstable_dict[conditions[11]]
    for k in range(1, 11):
        total_sg_genes = union(total_sg_genes, sg_cv_unstable_dict[conditions[k]])
        total_mg_genes = union(total_mg_genes, mg_cv_unstable_dict[conditions[k + 11]])

    print("SG total unstable genes (CV = %f): %i" % (n_filter, len(total_sg_genes)))
    print("MG total unstable genes (CV = %f): %i" % (n_filter, len(total_mg_genes)))

    ##### Randomization
    if randomization:
        for j in range(0,100000):
            n = random.randint(1,len(em.lines)-1)
            line = em.lines[n].replace("\n","").split("\t")
            evalues = line[1::]
            evalues = [float(x) for x in evalues]
            sample = random.sample(evalues, 5)
            if any([evalue >= 5 for evalue in sample]):
                cv = variation(sample)
                cv_list_sg.append(cv)
                condition_list_sg.append("Randomization")
                exposure_list_sg.append("random")

        for j in range(0,100000):
            n = random.randint(1,len(em.lines)-1)
            line = em.lines[n].replace("\n","").split("\t")
            evalues = line[1::]
            evalues = [float(x) for x in evalues]
            sample = random.sample(evalues,3)
            if any([evalue >= 5 for evalue in sample]):
                cv = variation(sample)
                cv_list_mg.append(cv)
                condition_list_mg.append("Randomization")
                exposure_list_mg.append("random")




    df_cv = pd.DataFrame.from_dict(cv_dict)
    df_cv.index = index
    df_cv_sg = df_cv[SG_conditions]
    df_cv_mg = df_cv[MG_conditions]


    #### ECDF plot

    # color palette
    colors = ["#136D03", "#0AA3EF", "#068ACB", "#0575AD", "#05608E", "#03486A", "#F00404", "#C20404", "#A50404", "#750303", "#4E0202", "#585958"]

    df_cv_ecdf_sg["Transcript"] = transcript_list_sg
    df_cv_ecdf_sg["Coefficient of Variation"] = cv_list_sg
    df_cv_ecdf_sg["Condition"] = condition_list_sg
    df_cv_ecdf_sg["Exposure"] = exposure_list_sg

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))

    sns.displot(df_cv_ecdf_sg, x="Coefficient of Variation", hue="Condition", kind="ecdf")
    plt.savefig("%s_SG_ECDF.png" %out_label)
    plt.close()

    df_cv_ecdf_mg["Transcript"] = transcript_list_mg
    df_cv_ecdf_mg["Coefficient of Variation"] = cv_list_mg
    df_cv_ecdf_mg["Condition"] = condition_list_mg
    df_cv_ecdf_mg["Exposure"] = exposure_list_mg

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))

    sns.displot(df_cv_ecdf_mg, x="Coefficient of Variation", hue="Condition", kind="ecdf")
    plt.savefig("%s_MG_ECDF.png" %out_label)
    plt.close()

    ### Violin Plots

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))
    plt.rcParams["figure.figsize"] = (10, 10)
    plt.rcParams['font.size'] = '18'
    sns.violinplot(data=df_cv_ecdf_sg, x="Coefficient of Variation", y="Condition")
    plt.tight_layout()
    plt.savefig("%s_SG_Violin.png" % out_label)
    plt.close()

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))
    plt.rcParams["figure.figsize"] = (10, 10)
    plt.rcParams['font.size'] = '18'
    sns.violinplot(data=df_cv_ecdf_mg, x="Coefficient of Variation", y="Condition")
    plt.tight_layout()
    plt.savefig("%s_MG_Violin.png" % out_label)
    plt.close()

    ##### KDE plot
    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))
    sns.displot(df_cv_ecdf_sg, x="Coefficient of Variation", hue="Condition", kind="kde", common_norm=False, common_grid=True,  bw_adjust=.1)
    plt.savefig("%s_SG_KDE_low_adjust.png" % out_label)
    plt.close()

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))
    sns.displot(df_cv_ecdf_mg, x="Coefficient of Variation", hue="Condition", kind="kde", common_norm=False, common_grid=True,  bw_adjust=.1)
    plt.savefig("%s_MG_KDE_low_adjust.png" % out_label)
    plt.close()

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))
    sns.displot(df_cv_ecdf_sg, x="Coefficient of Variation", hue="Condition", kind="kde", common_norm=False,
                common_grid=True)
    plt.savefig("%s_SG_KDE.png" % out_label)
    plt.close()

    sns.set_theme()
    sns.set_palette(sns.color_palette(colors))
    sns.displot(df_cv_ecdf_mg, x="Coefficient of Variation", hue="Condition", kind="kde", common_norm=False,
                common_grid=True)
    plt.savefig("%s_MG_KDE.png" % out_label)
    plt.close()

    ## Expression vs CV plot

    df_cv_ecdf_sg["Average FPKM"] = expression_list_sg
    df_cv_ecdf_mg["Average FPKM"] = expression_list_mg

    df_cv_ecdf_sg["Average log(FPKM)"] = [log10(x) for x in expression_list_sg]
    df_cv_ecdf_mg["Average log(FPKM)"] = [log10(x) for x in expression_list_mg]

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_sg, x="Coefficient of Variation", y="Average FPKM",
                      kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_SG_EXvsCV.png" % out_label)
    plt.close()

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_mg, x="Coefficient of Variation", y="Average FPKM",
                      kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_MG_EXvsCV.png" % out_label)
    plt.close()

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_sg, x="Coefficient of Variation", y="Average log(FPKM)",
                    kind="kde", color="r", figsize=(5,10))
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_SG_EXlogvsCV.png" % out_label)
    plt.close()

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_mg, x="Coefficient of Variation", y="Average log(FPKM)",
                    kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)

    plt.savefig("%s_MG_EXlogvsCV.png" % out_label)
    print("hey")
    plt.close()

    ###  Expression vs CV plot desglosed by first and second

    df_cv_ecdf_sg_first = df_cv_ecdf_sg[df_cv_ecdf_sg["Exposure"] == "first"]
    df_cv_ecdf_sg_second = df_cv_ecdf_sg[df_cv_ecdf_sg["Exposure"] == "second"]
    df_cv_ecdf_mg_first = df_cv_ecdf_mg[df_cv_ecdf_mg["Exposure"] == "first"]
    df_cv_ecdf_mg_second = df_cv_ecdf_mg[df_cv_ecdf_mg["Exposure"] == "second"]
    print(df_cv_ecdf_sg_first)

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_sg_first, x="Coefficient of Variation", y="Average log(FPKM)",
                      kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_SG_EXlogvsCV_first.png" % out_label)
    plt.close()

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_mg_first, x="Coefficient of Variation", y="Average log(FPKM)",
                      kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_MG_EXlogvsCV_first.png" % out_label)
    plt.close()

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_sg_second, x="Coefficient of Variation", y="Average log(FPKM)",
                      kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_SG_EXlogvsCV_second.png" % out_label)
    plt.close()

    sns.set_theme()
    g = sns.jointplot(data=df_cv_ecdf_mg_second, x="Coefficient of Variation", y="Average log(FPKM)",
                      kind="kde", color="r")
    g.plot_joint(sns.scatterplot, s=5)
    plt.savefig("%s_MG_EXlogvsCV_second.png" % out_label)
    plt.close()

    #
    #
    # ### Get 200 more expressed and unstable, histogram of conditions
    #
    # df_cv_ecdf_sg_unstable = df_cv_ecdf_sg[df_cv_ecdf_sg["Coefficient of Variation"] > 1.9]
    # df_cv_ecdf_mg_unstable = df_cv_ecdf_mg[df_cv_ecdf_sg["Coefficient of Variation"] > 1.35]
    #
    # df_cv_ecdf_sg_unstable_byexp = df_cv_ecdf_sg_unstable.sort_values("Average log(FPKM)", ascending=False)
    # df_cv_ecdf_mg_unstable_byexp = df_cv_ecdf_mg_unstable.sort_values("Average log(FPKM)", ascending=False)
    #
    # df_cv_ecdf_sg_unstable_byexp_200 = df_cv_ecdf_sg_unstable_byexp.head(200)
    # df_cv_ecdf_mg_unstable_byexp_200 = df_cv_ecdf_mg_unstable_byexp.head(200)
    #
    # df_cv_ecdf_sg_unstable_byexp_200.to_csv("%s_unstable_byexp_SG.tsv" %out_label, sep="\t")
    # df_cv_ecdf_mg_unstable_byexp_200.to_csv("%s_unstable_byexp_MG.tsv" % out_label, sep="\t")
    #
    # sns.set_theme()
    # sns.countplot(data=df_cv_ecdf_sg_unstable_byexp_200, x="Condition", order=SG_conditions, color="b")
    # plt.xticks(rotation=90)
    # plt.tight_layout()
    # plt.savefig("%s_unstable_byexp_SG_conditionhist.png" % out_label)
    # plt.close()
    #
    # sns.set_theme()
    # sns.set_palette(sns.color_palette(colors))
    # sns.countplot(data=df_cv_ecdf_mg_unstable_byexp_200, x="Condition", order=MG_conditions, color="b")
    # plt.xticks(rotation=90)
    # plt.tight_layout()
    # plt.savefig("%s_unstable_byexp_MG_conditionhist.png" % out_label)
    # plt.close()
    #
    #
    #
    #
    #
    # ######
    #
    #
    # ax = df_cv.plot.kde()
    # ax.set_xlabel("Coefficient of variation")
    # fig = ax.get_figure()
    # fig.savefig("%s_distribution.pdf" % out_label)
    #
    # for condition in conditions:
    #     df_cv["%s_rank" %condition] = df_cv[condition].rank(method="first")
    #     rank_df["%s_rank" %condition] = df_cv[condition].rank(method="first")
    #     if condition.startswith("MG"):
    #         mg_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #         if condition.startswith("MG unfed"):
    #             mg_unfed_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #         elif condition.startswith("MG first"):
    #             mg_1_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #         elif condition.startswith("MG second"):
    #             mg_2_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #     elif condition.startswith("SG"):
    #         sg_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #         if condition.startswith("SG unfed"):
    #             sg_unfed_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #         elif condition.startswith("SG first"):
    #             sg_1_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    #         elif condition.startswith("SG second"):
    #             sg_2_rank_df["%s_rank" % condition] = df_cv[condition].rank(method="first")
    # df_cv["Rank_mean"] = rank_df.mean(axis=1, skipna=True)
    # df_cv["SG_rank_mean"] = sg_rank_df.mean(axis=1, skipna=True)
    # df_cv["MG_rank_mean"] = mg_rank_df.mean(axis=1, skipna=True)
    # df_cv["SG_Unfed_rank_mean"] = sg_unfed_rank_df.mean(axis=1, skipna=True)
    # df_cv["SG_first_rank_mean"] = sg_1_rank_df.mean(axis=1, skipna=True)
    # df_cv["SG_second_rank_mean"] = sg_2_rank_df.mean(axis=1, skipna=True)
    # df_cv["MG_unfed_rank_mean"] = mg_unfed_rank_df.mean(axis=1, skipna=True)
    # df_cv["MG_first_rank_mean"] = mg_1_rank_df.mean(axis=1, skipna=True)
    # df_cv["MG_second_rank_mean"] = mg_2_rank_df.mean(axis=1, skipna=True)
    # df_cv.to_csv(output, sep="\t")
    #
    # ax = df_cv_sg.plot.kde()
    # ax.set_xlabel("Coefficient of variation")
    # fig = ax.get_figure()
    # fig.savefig("%s_sg_distribution.pdf" % out_label)
    # ax = df_cv_mg.plot.kde()
    # ax.set_xlabel("Coefficient of variation")
    # fig = ax.get_figure()
    # fig.savefig("%s_mg_distribution.pdf" % out_label)
    #
    # ## File for https://network.shinyapps.io/superexacttest/
    # df_cv_unstable_sg = pd.DataFrame.from_dict(sg_cv_unstable_dict, orient="index").T
    # df_cv_unstable_mg = pd.DataFrame.from_dict(mg_cv_unstable_dict, orient="index").T
    #
    # df_cv_unstable_sg.to_csv("%s_unstable_cv%f_SG.tsv" %(out_label,n_filter), sep="\t", index=False)
    # df_cv_unstable_mg.to_csv("%s_unstable_cv%f_MG.tsv" %(out_label,n_filter), sep="\t", index=False)
    #
    # # with open("%s_unstable_SG.tsv" %out_label, "wt") as out:
    # #     for condition in sg_cv_unstable_dict:
    # #         out.write("%s\t%s\n" %(condition,"\t".join(sg_cv_unstable_dict[condition])))
    # # with open("%s_unstable_MG.tsv" %out_label, "wt") as out:
    # #     for condition in mg_cv_unstable_dict:
    # #         out.write("%s\t%s\n" %(condition,"\t".join(mg_cv_unstable_dict[condition])))
    #
    # ###1000 more unstable.
    #
    # sg_cv_unstable_1000_dict = dict()
    # mg_cv_unstable_1000_dict = dict()
    # for k in range(0,11):
    #     sg_cv_unstable_1000_dict[conditions[k]] = []
    #     mg_cv_unstable_1000_dict[conditions[k+11]] = []
    #     cvs_sg = [x for x in cv_tuple_dict[conditions[k]] if str(lambda x: x[1]) != 'nan']
    #     cvs_sg.sort(key=lambda x: x[1], reverse=True)
    #     cvs_mg = [x for x in cv_tuple_dict[conditions[k+11]] if str(lambda x: x[1]) != 'nan']
    #     cvs_mg.sort(key=lambda x: x[1], reverse=True)
    #     for j in range(0,m_filter):
    #         sg_cv_unstable_1000_dict[conditions[k]].append(cvs_sg[j][0])
    #         mg_cv_unstable_1000_dict[conditions[k+11]].append(cvs_mg[j][0])
    #
    #
    #
    # df_cv_unstable_sg = pd.DataFrame.from_dict(sg_cv_unstable_1000_dict, orient="index").T
    # df_cv_unstable_mg = pd.DataFrame.from_dict(mg_cv_unstable_1000_dict, orient="index").T
    #
    # df_cv_unstable_sg.to_csv("%s_unstable_%i_SG.tsv" % (out_label, m_filter), sep="\t", index=False)
    # df_cv_unstable_mg.to_csv("%s_unstable_%i_MG.tsv" % (out_label, m_filter), sep="\t", index=False)
    #
    # #calculation of the background
    # total_sg_genes = sg_cv_unstable_1000_dict[conditions[0]]
    # total_mg_genes = mg_cv_unstable_1000_dict[conditions[11]]
    # for k in range(1,11):
    #     total_sg_genes = union(total_sg_genes, sg_cv_unstable_1000_dict[conditions[k]])
    #     total_mg_genes = union(total_mg_genes, mg_cv_unstable_1000_dict[conditions[k+11]])
    #
    # print("SG total unstable genes (%i per sample): %i" %(m_filter, len(total_sg_genes)))
    # print("MG total unstable genes (%i per sample): %i" %(m_filter, len(total_mg_genes)))
    #
    # ###1000 more stable.
    #
    # sg_cv_stable_1000_dict = dict()
    # mg_cv_stable_1000_dict = dict()
    # for k in range(0, 11):
    #     sg_cv_stable_1000_dict[conditions[k]] = []
    #     mg_cv_stable_1000_dict[conditions[k + 11]] = []
    #     cvs_sg = [x for x in cv_tuple_dict[conditions[k]] if str(lambda x: x[1]) != 'nan']
    #     cvs_sg.sort(key=lambda x: x[1])
    #     cvs_mg = [x for x in cv_tuple_dict[conditions[k + 11]] if str(lambda x: x[1]) != 'nan']
    #     cvs_mg.sort(key=lambda x: x[1])
    #     for j in range(0, m_filter):
    #         sg_cv_stable_1000_dict[conditions[k]].append(cvs_sg[j][0])
    #         mg_cv_stable_1000_dict[conditions[k + 11]].append(cvs_mg[j][0])
    #
    # df_cv_stable_sg = pd.DataFrame.from_dict(sg_cv_stable_1000_dict, orient="index").T
    # df_cv_stable_mg = pd.DataFrame.from_dict(mg_cv_stable_1000_dict, orient="index").T
    #
    # df_cv_stable_sg.to_csv("%s_stable_%i_SG.tsv" % (out_label, m_filter), sep="\t", index=False)
    # df_cv_stable_mg.to_csv("%s_stable_%i_MG.tsv" % (out_label, m_filter), sep="\t", index=False)
    #
    # # calculation of the background
    # total_sg_genes = sg_cv_stable_1000_dict[conditions[0]]
    # total_mg_genes = mg_cv_stable_1000_dict[conditions[11]]
    # for k in range(1, 11):
    #     total_sg_genes = union(total_sg_genes, sg_cv_stable_1000_dict[conditions[k]])
    #     total_mg_genes = union(total_mg_genes, mg_cv_stable_1000_dict[conditions[k + 11]])
    #
    # print("SG total stable genes (%i per sample): %i" % (m_filter, len(total_sg_genes)))
    # print("MG total stable genes (%i per sample): %i" % (m_filter, len(total_mg_genes)))

def union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list

def intersection(lst1,lst2,reverse=False):
    if reverse:
        lst3 = [value for value in lst1 if value not in lst2]
    else:
        lst3 = [value for value in lst1 if value in lst2]
    return lst3

def freq(lst):
    d = {}
    for i in lst:
        if d.get(i):
            d[i] += 1
        else:
            d[i] = 1
    return d

def order_dictionary(D):
    D_ord = {}
    key_ord = sorted(D, key=D.get)
    for key in reversed(key_ord):
        D_ord[key] = D[key]
    return D_ord

def filter_dictionary(D, n):
    D_filt = {}
    for key in D:
        if float(D[key]) >= n:
            D_filt[key] = D[key]
    return D_filt

def piechart(labels, sizes, out):
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    # labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    # sizes = [15, 30, 45, 10]
    # explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    sns.set_theme(style="whitegrid")

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    fig1.savefig(out)

def barplot_from_dict(D,out,title, label):
    values = list(D.values())
    keys = list(D.keys())
    keys = [x.replace("\n", "") for x in keys]
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 6))
    plt.title(title, fontsize=20)
    plt.bar(range(len(D)), values, align='center')
    plt.xticks(range(len(D)), keys, rotation='vertical')
    plt.ylabel(label)
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.3)
    plt.savefig(out)
    plt.close()

def barplot_from_dict_error(D, D_error, out,title, label):
    values = list(D.values())
    keys = list(D.keys())
    keys = [x.replace("\n", "") for x in keys]
    error = list(D_error.values())
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 6))
    plt.title(title, fontsize=20)
    plt.bar(range(len(D)), values, align='center')
    plt.xticks(range(len(D)), keys, rotation='vertical')
    plt.errorbar(range(len(D)), values, error, fmt='.', color='Black', elinewidth=2, capthick=5, errorevery=1, alpha=0.5, ms=4,
                 capsize=2)
    plt.ylabel(label)
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.3)
    plt.savefig(out, dpi=300)
    plt.close()

def barhplot_from_dict(D,out,label,number, filter_n=False):
    values = list(D.values())
    values.reverse()
    keys = list(D.keys())
    keys.reverse()
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(15, 10))
    plt.barh(range(len(D)), values, align='center', fontsize=20)
    plt.yticks(range(len(D)), keys)
    plt.xlabel(label)
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.3)
    plt.savefig(out, dpi=300)
    plt.close()
    if len(D) > number:
        values = list(D.values())[0:number]
        values.reverse()
        keys = list(D.keys())[0:number]
        keys.reverse()
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(15, 10))
        plt.barh(range(number), values, align='center', fontsize=20)
        plt.yticks(range(number), keys)
        plt.xlabel(label)
        plt.tight_layout()
        # plt.subplots_adjust(bottom=0.3)
        plt.savefig("%s_%i.tiff" %(out.replace(".tiff",""), number), dpi=300)
        plt.close()
    else:
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(15, 10))
        plt.barh(range(len(D)), values, align='center', fontsize=20)
        plt.yticks(range(len(D)), keys)
        plt.xlabel(label)
        plt.tight_layout()
        # plt.subplots_adjust(bottom=0.3)
        plt.savefig("%s_%i.tiff" %(out.replace(".tiff",""), number), dpi=300)
        plt.close()

def weight_figure(output_label):
    rename_dict = dict()
    rename_dict["pre_2_SG_S235"] = "SG_unfed_1,5mg_pre2"
    rename_dict["pre_4_SG_S267"] = "SG_unfed_1,7mg_pre4"
    rename_dict["pre_6_SG_S236"] = "SG_unfed_1,4mg_pre6"
    rename_dict["pre_8_SG_S274"] = "SG_unfed_1,6mg_pre8"
    rename_dict["pre_10_SG_S237"] = "SG_unfed_1,4mg_pre10"
    rename_dict["1_12h_SG_S268"] = "SG_1_12h_1,5mg_1"
    rename_dict["105_12h_SG_S1"] = "SG_1_12h_2mg_105"
    rename_dict["5_12h_SG_S224"] = "SG_1_12h_1,8mg_5"
    rename_dict["101_12h_SG_S225"] = "SG_1_12h_1,7mg_101"
    rename_dict["103_12h_SG_S275"] = "SG_1_12h_2mg_103"
    rename_dict["7_24h_SG_S270"] = "SG_1_24h_1,5mg_7"
    rename_dict["9_24h_SG_S226"] = "SG_1_24h_1,7mg_9"
    rename_dict["109_24h_SG_S312"] = "SG_1_24h_2,3mg_109"
    rename_dict["106_24h_SG_S227"] = "SG_1_24h_1,6mg_106"
    rename_dict["8_24h_SG_S269"] = "SG_1_24h_2,3mg_8"
    rename_dict["13_48h_SG_S255"] = "SG_1_48h_3,2mg_13"
    rename_dict["15_48h_SG_S256"] = "SG_1_48h_3mg_15"
    rename_dict["111_48h_SG_S228"] = "SG_1_48h_2,9mg_111"
    rename_dict["113_48h_SG_S229"] = "SG_1_48h_3,1mg_113"
    rename_dict["114_48h_SG_S271"] = "SG_1_48h_2,9mg_114"
    rename_dict["18_72h_SG_S257"] = "SG_1_72h_5,2mg_18"
    rename_dict["19_72h_SG_S258"] = "SG_1_72h_7,4mg_19"
    rename_dict["20_72h_SG_S259"] = "SG_1_72h_5,5mg_20"
    rename_dict["116_72h_SG_S307"] = "SG_1_72h_6,2mg_116"
    rename_dict["118_72h_SG_S230"] = "SG_1_72h_5,9mg_118"
    rename_dict["21_96h_SG_S239"] = "SG_1_96h_10,3mg_21"
    rename_dict["22_96h_SG_S240"] = "SG_1_96h_10,3mg_22"
    rename_dict["24_96h_SG_S260"] = "SG_1_96h_10,2mg_24"
    rename_dict["121_96h_SG_S231"] = "SG_1_96h_9,1mg_121"
    rename_dict["122_96h_SG_S261"] = "SG_1_96h_9,5mg_122"
    rename_dict["52_12h_SG_S280"] = "SG_2_12h_2,2mg_52"
    rename_dict["53_12h_SG_S232"] = "SG_2_12h_1,9mg_53"
    rename_dict["152_12h_SG_S233"] = "SG_2_12h_2mg_152"
    rename_dict["153_12h_SG_S234"] = "SG_2_12h_1,7mg_153"
    rename_dict["155_12h_SG_S272"] = "SG_2_12h_1,8mg_155"
    rename_dict["56_1d_SG_S262"] = "SG_2_24h_2,2mg_56"
    rename_dict["59_1d_SG_S263"] = "SG_2_24h_2,4mg_59"
    rename_dict["156_1d_SG_S241"] = "SG_2_24h_2,2mg_156"
    rename_dict["157_1d_SG_S273"] = "SG_2_24h_2,3mg_157"
    rename_dict["159_1d_SG_S242"] = "SG_2_24h_2,5mg_159"
    rename_dict["62_2d_SG_S243"] = "SG_2_48h_4,6mg_62"
    rename_dict["63_2d_SG_S244"] = "SG_2_48h_4,5mg_63"
    rename_dict["65_2d_SG_S245"] = "SG_2_48h_4,2mg_65"
    rename_dict["162_2d_SG_S246"] = "SG_2_48h_5mg_162"
    rename_dict["163_2d_SG_S264"] = "SG_2_48h_4,3mg_163"
    rename_dict["67_3d_SG_S247"] = "SG_2_72h_8,7mg_67"
    rename_dict["68_3d_SG_S265"] = "SG_2_72h_8,2mg_68"
    rename_dict["69_3d_SG_S248"] = "SG_2_72h_9,3mg_69"
    rename_dict["166_3d_SG_S249"] = "SG_2_72h_6,7mg_166"
    rename_dict["167_3d_SG_S250"] = "SG_2_72h_7,2mg_167"
    rename_dict["71_4d_SG_S266"] = "SG_2_96h_11,1mg_71"
    rename_dict["74_4d_SG_S251"] = "SG_2_96h_12,9mg_74"
    rename_dict["75_4d_SG_S252"] = "SG_2_96h_15,1mg_75"
    rename_dict["171_4d_SG_S253"] = "SG_2_96h_10,2mg_171"
    rename_dict["172_4d_SG_S254"] = "SG_2_96h_10,1mg_172"
    rename_dict["pre4_MG_S24"] = "MG_unfed_1,7mg_pre4"
    rename_dict["pre10_MG_S2"] = "MG_unfed_1,4mg_pre10"
    rename_dict["pre_8_MG_S279"] = "MG_unfed_1,6mg_pre8"
    rename_dict["4_12h_MG_S310"] = "MG_1_12h_1,8mg_4"
    rename_dict["5_12h_MG_S311"] = "MG_1_12h_1,8mg_5"
    rename_dict["104_12h_MG_S238"] = "MG_1_12h_1,8mg_104"
    rename_dict["9_24h_MG_S297"] = "MG_1_24h_1,7mg_9"
    rename_dict["10_24h_MG_S298"] = "MG_1_24h_1,6mg_10"
    rename_dict["107_24h_MG_S299"] = "MG_1_24h_1,9mg_107"
    rename_dict["13_48h_MG_S281"] = "MG_1_48h_3,2mg_13"
    rename_dict["15_48h_MG_S282"] = "MG_1_48h_3mg_15"
    rename_dict["114_48h_MG_S283"] = "MG_1_48h_2,9mg_114"
    rename_dict["20_72h_MG_S300"] = "MG_1_72h_5,5mg_20"
    rename_dict["116_72h_MG_S284"] = "MG_1_72h_6,2mg_116"
    rename_dict["118_72h_MG_S285"] = "MG_1_72h_5,9mg_118"
    rename_dict["22_96h_MG_S286"] = "MG_1_96h_10,3mg_22"
    rename_dict["24_96h_MG_S287"] = "MG_1_96h_10,2mg_24"
    rename_dict["122_96h_MG_S288"] = "MG_1_96h_9,5mg_122"
    rename_dict["155_12h_MG_S308"] = "MG_2_12h_1,8mg_155"
    rename_dict["54_12_MG_S309"] = "MG_2_12h_2,2mg_54"
    rename_dict["152_12h_MG_S276"] = "MG_2_12h_2mg_152"
    rename_dict["56_1d_MG_S277"] = "MG_2_24h_2,2mg_56"
    rename_dict["59_1d_MG_S301"] = "MG_2_24h_2,4mg_59"
    rename_dict["157_1d_MG_S278"] = "MG_2_24h_2,3mg_157"
    rename_dict["62_2d_MG_S289"] = "MG_2_48h_4,6mg_62"
    rename_dict["63_2d_MG_S290"] = "MG_2_48h_4,5mg_63"
    rename_dict["163_2d_MG_S302"] = "MG_2_48h_4,3mg_163"
    rename_dict["67_3d_MG_S291"] = "MG_2_72h_8,7mg_67"
    rename_dict["68_3d_MG_S292"] = "MG_2_72h_8,2mg_68"
    rename_dict["167_3d_MG_S293"] = "MG_2_72h_7,2mg_167"
    rename_dict["71_4d_MG_S294"] = "MG_2_96h_11,1mg_71"
    rename_dict["74_4d_MG_S295"] = "MG_2_96h_12,9mg_74"
    rename_dict["171_4d_MG_S296"] = "MG_2_96h_10,2mg_171"

    weights = []
    conditions = []
    exposures = []
    hours = []

    for key in rename_dict:
        tissue = rename_dict[key].split("_")[0]
        if tissue == "SG":
            weights.append(float(rename_dict[key].split("mg_")[0].split("_")[-1].replace(",",".")))
            exp =  rename_dict[key].split("_")[1]
            if exp == "1":
                exp = "first exposure"
            elif exp == "2":
                exp="second exposure"
            if exp == "unfed":
                hour = "0h"
                condition = "A unfed"
            else:
                hour = rename_dict[key].split("_")[2]
                condition = exp + " " + hour
            conditions.append(condition)
            exposures.append(exp)
            hours.append(hour)

    df = pd.DataFrame(data={"Weight (mg)":weights, "Condition":conditions, "Exposure":exposures, "Feeding time":hours})
    df = df.sort_values(by=['Condition'])

    sns.set_theme()
    ax = sns.catplot(x="Condition", y="Weight (mg)", data=df)
    plt.xticks(rotation=90)
    ax.savefig("%s_weightdistribution.png" %output_label)
    plt.close()

    df_dynamic = df[df["Exposure"] != "unfed"]
    sns.set_theme()
    ax =sns.catplot(x="Feeding time", y="Weight (mg)", hue="Exposure", kind="point", data=df_dynamic)
    plt.xticks(rotation=40, ha="right")
    ax.savefig("%s_weightdynamic.png" %output_label, dpi=300)
    plt.close()

def exosome_ematrix_analysis(ematrix, out_label, full=False):

    em = Ematrix(ematrix)
    header = em.lines[0].replace("\n","").split("\t")
    color_vector = []
    form_vector = []
    for k in range(1, len(header)):
        sample = header[k]
        if "FTS" in sample:
            color_vector.append(header[k].split("-")[1].split("_")[0])
            form_vector.append("FTS")
        elif "Undetermined" in sample:
            color_vector.append("Undetermined")
            form_vector.append("Undetermined")
        elif "Michael" in sample:
            color_vector.append("Bilbao")
            form_vector.append("Bilbao")
        else:
            if sample.startswith("1"):
                color_vector.append(header[k].split("-")[1].split("_")[0])
                form_vector.append("First exposure")
            elif sample.startswith("2"):
                color_vector.append(header[k].split("-")[1].split("_")[0])
                form_vector.append("Second exposure")

    em.PCA(out_label, color_vector, form_vector)
    em.MDS(out_label, color_vector, form_vector)
    em.correlation_matrix(out_label)
    em.hierarchical_clustering(out_label, "all")

def exosome_alldata_ematrix_analysis(ematrix, samplesheet, out_label):
    em = Ematrix(ematrix)
    ss = SampleSheet(samplesheet)

    header = getrow(em.lines[0])
    color_vector = []
    form_vector = []

    exp1_2019 = ["L1","L3","L5","L7","L10","L12"]
    exp2_2019 = ["L2","L4","L6","L8", "L11", "L13"]


    for k in range(1, len(header)):
        sample = header[k]

        if ss.dict[sample] == "2017":
            color_vector.append("2017")
            if "MG" in sample:
                form_vector.append("MG")
            elif "SG" in sample:
                form_vector.append("SG")
            else:
                form_vector.append("Other")
        elif ss.dict[sample] == "2019":
            color_vector.append("2019")
            n = sample.split("_")[1]
            if n in exp1_2019:
                form_vector.append("1EXP")
            elif n in exp2_2019:
                form_vector.append("2EXP")
            else:
                form_vector.append("Other")
        elif ss.dict[sample] == "2023":
            color_vector.append("2023")
            if sample.startswith("1"):
                form_vector.append("1EXP")
            elif sample.startswith("2"):
                form_vector.append("2EXP")
            else:
                form_vector.append("Other")

    em.PCA(out_label, color_vector, form_vector)
    em.MDS(out_label, color_vector, form_vector)
    em.correlation_matrix(out_label)
    em.hierarchical_clustering(out_label, "all")

def virus_segmentation(file, output):
    fa = Fasta(file)
    with open(output, "wt") as out:
        for key in fa.dict:
            parts = math.ceil(len(fa.dict[key]) / 1000)
            bases = round(len(fa.dict[key]) / parts)

            for k in range(0,parts):
                if k != parts - 1:
                    seq = fa.dict[key][k * bases : (k+1) * bases]
                else:
                    seq = fa.dict[key][k * bases ::]
                out.write(">%s_%i\n%s\n" % (key, k + 1, seq))

def ncRNA_stats(gff, output):
    gff = Gff(gff)
    stats = {"Transfer RNA":0,"Ribosomal RNA":0, "Signal recognition particle":0, "Ribozyme":0, "Small nuclear RNA":0,
                       "Small nucleolar RNA":0, "Binding small RNA":0, "MicroRNA":0, "Long non-coding RNA":10696, "UTR stem-loop":0,
                       "Riboswitch":0, "Iron response element":0, "Potassium channel RNA editing signal":0}
    for line in gff.lines:
        gff_line = Gff.Line(line)
        if gff_line.type=="gene":
            if "ribosomal_rna" in line.lower():
                stats["Ribosomal RNA"] += 1
            elif "signal_recognition_particle" in line.lower():
                stats["Signal recognition particle"] += 1
            elif "trna" in line.lower():
                stats["Transfer RNA"] += 1
            elif "rnase" in line.lower():
                stats["Ribozyme"] += 1
            elif "spliceosomal_rna" in line.lower():
                stats["Small nuclear RNA"] += 1
            elif "small_nucleolar" in line.lower():
                stats["Small nucleolar RNA"] += 1
            elif "binding_rna" in line.lower():
                stats["Binding small RNA"] += 1
            elif "mirna" in line.lower():
                stats["MicroRNA"] += 1
            elif "utr_stem-loop" in line.lower():
                stats["UTR stem-loop"] += 1
            elif "riboswitch" in line.lower():
                stats["Riboswitch"] += 1
            elif "iron_response" in line.lower():
                stats["Iron response element"] += 1
            elif "potassium" in line.lower():
                stats["Potassium channel RNA editing signal"] += 1

    with open(output, "wt") as out:
        out.write("Type\tClass\tNumber\n"
                  "Transfer RNA\tNon-coding RNA\t%i\n"
                  "Ribosomal RNA\tNon-coding RNA\t%i\n"
                  "Small nuclear RNA\tNon-coding RNA\t%i\n"
                  "Small nucleolar RNA\tNon-coding RNA\t%i\n"
                  "Binding small RNA\tNon-coding RNA\t%i\n"
                  "microRNA\tNon-coding RNA\t%i\n"
                  "Long non-coding RNA\tNon-coding RNA\t%i\n"
                  "Signal recognition particle\tNon-coding RNA\t%i\n"
                  "Ribozyme\tNon-coding RNA\t%i\n"
                  "Riboswitch\tCis-regulatory element\t%i\n"
                  "UTR stem-loop\tCis-regulatory element\t%i\n"
                  "Iron response element\tCis-regulatory element\t%i\n"
                  "Potassium chanel RNA editing signal\tCis-regulatory element\t%i\n"
                  %(stats["Transfer RNA"], stats["Ribosomal RNA"], stats["Small nuclear RNA"], stats["Small nucleolar RNA"], stats["Binding small RNA"],
                    stats["MicroRNA"], stats["Long non-coding RNA"], stats["Signal recognition particle"], stats["Ribozyme"],
                    stats["Riboswitch"], stats["UTR stem-loop"], stats["Iron response element"], stats["Potassium channel RNA editing signal"])
                  )



def process_line(line, ids):
    fields = line.strip().split('\t')
    attributes = fields[8].split(';')

    id_field = attributes[0].replace("ID=","")
    if 'exon' in id_field:
        exon_id = ".".join(id_field.split(".")[0:2])
    else:
        exon_id = False

    if exon_id and exon_id in ids and float(fields[5]) < ids[exon_id]:
        attributes = [x + ".2" for x in attributes]
    else:
        ids[exon_id] = float(fields[5])


    if id_field in ids and float(fields[5]) < ids[id_field]:
        attributes = [x + ".2" for x in attributes]

    else:
        ids[id_field] = float(fields[5])

    fields[8] = ';'.join(attributes)
    return '\t'.join(fields)



def get_MGSG_means(ematrix, output):
    df = pd.read_csv(ematrix, sep='\t', skiprows=[1])
    # Group columns based on SG or MG
    sg_columns = [col for col in df.columns if 'SG' in col]
    mg_columns = [col for col in df.columns if 'MG' in col]

    # Calculate mean expression for SG and MG
    df['Mean_Expression_SG'] = df[sg_columns].mean(axis=1)
    df['Mean_Expression_MG'] = df[mg_columns].mean(axis=1)

    df.to_csv(output, sep="\t")










