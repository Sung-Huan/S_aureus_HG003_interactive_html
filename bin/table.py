from datetime import date
from random import randint
from gff3 import Gff3Parser
from bokeh.io import output_file, save
from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models import (ColumnDataSource, DataRange1d, Plot, LinearAxis,
                          Grid, Circle, HoverTool, BoxSelectTool, LabelSet,
                          Label)
from bokeh.models.widgets import (DataTable, TableColumn, DateFormatter,
                                  StringFormatter, NumberFormatter,
                                  HTMLTemplateFormatter, StringEditor,
                                  IntEditor, NumberEditor, SelectEditor)
from bokeh.sampledata.periodic_table import elements
from bokeh.plotting import figure
from bokeh.resources import INLINE
import argparse
import math
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="gene annotation gff file")
parser.add_argument("-n","--name_file",help="gene name file")
parser.add_argument("-g","--gene_quanti_file",help="gene quantification file")
parser.add_argument("-f","--features",help="input features")
args = parser.parse_args()

def check_ps(cds, gffs, type_, data):
    detect = False
    tss_list = []
    for tss_ps in gffs:
        if cds.strand == tss_ps.strand:
            if (cds.start <= (tss_ps.start + 3)) and (
                    cds.end >= (tss_ps.start - 3)):
                detect = True
                tss_list.append("-".join([str(tss_ps.start),
                                          str(tss_ps.end)]))
    if detect:
        data[type_].append("\n".join(tss_list))
    else:
        data[type_].append("NA")

def check_tss(entry, gffs, type_, data, feature):
    detect = False
    tss_list = []
    for tss_ps in gffs:
        if feature == "CDS":
            if "associated_gene" in tss_ps.attributes.keys():
                if (entry.attributes["locus_tag"] in 
                    tss_ps.attributes["associated_gene"]):
                    detect = True
                    tss_list.append("-".join([str(tss_ps.start),
                                              str(tss_ps.end)]))
        else:
            if entry.strand == tss_ps.strand:
                if (entry.start <= (tss_ps.start + 3)) and (
                        entry.end >= (tss_ps.start - 3)):
                    detect = True
                    tss_list.append("-".join([str(tss_ps.start),
                                              str(tss_ps.end)]))
            else:
                if ((entry.start - 100) <= tss_ps.start) and (
                        (entry.end + 100) >= tss_ps.start):
                    detect = True
                    tss_list.append("-".join([str(tss_ps.start),
                                              str(tss_ps.end)]))
    if detect:
        data[type_].append("\n".join(tss_list))
    else:
        data[type_].append("NA")

def check_term(cds, terms, in_f, out_f, data):
    detect = False
    term_list = []
    for term in terms:
        if term.strand == cds.strand:
            if term.strand == "+":
                start = cds.end - in_f
                if start < cds.start:
                    start = cds.start
                end = cds.end + out_f
            else:
                start = cds.start - out_f
                end = cds.start + in_f
                if end > cds.end:
                    end = cds.end
            if ((term.start >= start) and (term.end <= end)) or (
                    (term.start <= start) and (term.end >= end)) or (
                    (term.start <= start) and (term.end >= start) and (
                     term.end <= end)) or (
                    (term.end >= end) and (term.start <= end) and (
                     term.start >= start)):
                term_list.append("-".join([str(term.start),
                                           str(term.end)]))
                detect = True
    if detect:
        data["term"].append("\n".join(term_list))
    else:
        data["term"].append("NA")

def get_quant(cds, quants, data, type_):
    for quant in quants:
        if quant[0] == "sense":
            if (type_ == quant[3]) and (
                    str(cds.start) == quant[4]) and (
                    str(cds.end) == quant[5]) and (
                    cds.strand == quant[7]):
                data["TSB_OD_0.2"].append(float(quant[10]))
                data["TSB_OD_0.2_TEX"].append(float(quant[11]))
                data["TSB_OD_0.5"].append(float(quant[12]))
                data["TSB_OD_0.5_TEX"].append(float(quant[13]))
                data["TSB_OD_1"].append(float(quant[14]))
                data["TSB_OD_1_TEX"].append(float(quant[15]))
                data["TSB_ON"].append(float(quant[16]))
                data["TSB_ON_TEX"].append(float(quant[17]))
                data["TSB_t0"].append(float(quant[18]))
                data["TSB_t0_TEX"].append(float(quant[19]))
                data["TSB_t1"].append(float(quant[20]))
                data["TSB_t1_TEX"].append(float(quant[21]))
                data["TSB_t2"].append(float(quant[22]))
                data["TSB_t2_TEX"].append(float(quant[23]))
                data["pMEM_OD_0.2"].append(float(quant[24]))
                data["pMEM_OD_0.2_TEX"].append(float(quant[25]))
                data["pMEM_OD_0.5"].append(float(quant[26]))
                data["pMEM_OD_0.5_TEX"].append(float(quant[27]))
                data["pMEM_OD_1"].append(float(quant[28]))
                data["pMEM_OD_1_TEX"].append(float(quant[29]))
                data["pMEM_ON"].append(float(quant[30]))
                data["pMEM_ON_TEX"].append(float(quant[31]))
                data["pMEM_t0"].append(float(quant[32]))
                data["pMEM_t0_TEX"].append(float(quant[33]))
                data["pMEM_t1"].append(float(quant[34]))
                data["pMEM_t1_TEX"].append(float(quant[35]))
                data["pMEM_t2"].append(float(quant[36]))
                data["pMEM_t2_TEX"].append(float(quant[37]))
                data["frag"].append(float(quant[38]))
                break

def check_tran(cds, trans, data):
    ta_list = []
    detect = False
    for tran in trans:
        if cds.strand == tran.strand:
            if ((cds.start >= tran.start) and (cds.end <= tran.end)) or (
                (cds.start <= tran.start) and (cds.end >= tran.end)) or (
                (cds.start >= tran.start) and (cds.start <= tran.end) and (
                 cds.end >= tran.end)) or (
                (cds.start <= tran.start) and (cds.end >= tran.start) and (
                 cds.end <= tran.end)):
                detect = True
                ta_list.append("-".join([str(tran.start), str(tran.end)]))
    if detect:
        data["tran"].append("\n".join(ta_list))
    else:
        data["tran"].append("NA")

def import_basic(data, entry, feature):
    data["features"].append(feature)
    data["start"].append(entry.start)
    data["end"].append(entry.end)
    data["strand"].append(entry.strand)

def get_cds_info(cdss, gffs, data, quants, names):
    for cds in cdss:
        detect = False
        import_basic(data, cds, "CDS")
        for locus, name in names.items():
            if cds.attributes["locus_tag"] == locus:
                if (name != "-") and (len(name) != 0):
                    data["gene_name"].append(name)
                else:
                    data["gene_name"].append(cds.attributes["locus_tag"])
        for gene in gffs["gene"]:
            if gene.attributes["ID"] in cds.attributes["Parent"].split(","):
                if "db_xref" in gene.attributes.keys():
                    link = ("https://www.ncbi.nlm.nih.gov/gene/" +
                            gene.attributes["db_xref"].split(":")[-1])
                    hyper = ('<a href= ' + link + ' target="_blank"> ' +
                             gene.attributes["db_xref"] + '</a>')
                    data["link"].append(hyper)
                else:
                    data["link"].append(gene.attributes["product"])
                detect = True
        if not detect:
            data["link"].append(cds.attributes["product"])
        check_tran(cds, gffs["transcript"], data)
        check_tss(cds, gffs["TSS"], "tss", data, "CDS")
        check_ps(cds, gffs["processing_site"], "ps", data)
        check_term(cds, gffs["terminator"], 30, 300, data)
        get_quant(cds, quants, data, "CDS")


def get_rfam_info(rfams, gffs, data, quants, feature):
    for rfam in rfams:
        import_basic(data, rfam, feature)
        if "rfam_id" in rfam.attributes.keys():
            link = ("http://rfam.xfam.org/family/" +
                    rfam.attributes["rfam_id"].split(":")[-1])
            hyper = ('<a href= ' + link + ' target="_blank"> ' +
                     rfam.attributes["rfam_id"] + '</a>')
            data["link"].append(hyper)
        else:
            data["link"].append(rfam.attributes["Name"])
        data["gene_name"].append(rfam.attributes["Name"])
        check_tran(rfam, gffs["transcript"], data)
        check_tss(rfam, gffs["TSS"], "tss", data, feature)
        check_ps(rfam, gffs["processing_site"], "ps", data)
        check_term(rfam, gffs["terminator"], 10, 30, data)
        get_quant(rfam, quants, data, feature)

def get_rna_info(rnas, gffs, data, quants, feature):
    for rna in rnas:
        detect = False
        import_basic(data, rna, feature)
        if "note" in rna.attributes.keys():
            data["gene_name"].append(rna.attributes["note"])
        else:
            data["gene_name"].append(rna.attributes["product"])
        if "db_xref" in rna.attributes.keys():
            link = ("https://www.ncbi.nlm.nih.gov/gene/" +
                    rna.attributes["db_xref"].split(":")[-1])
            hyper = ('<a href= ' + link + ' target="_blank"> ' +
                     rna.attributes["db_xref"] + '</a>')
            data["link"].append(hyper)
        else:
            data["link"].append(rna.attributes["product"])
        check_tran(rna, gffs["transcript"], data)
        check_tss(rna, gffs["TSS"], "tss", data, feature)
        check_ps(rna, gffs["processing_site"], "ps", data)
        check_term(rna, gffs["terminator"], 30, 30, data)
        get_quant(rna, quants, data, feature)

def get_crispr_info(crisprs, gffs, data, quants):
    for crispr in crisprs:
        import_basic(data, crispr, "CRISPR")
        data["gene_name"].append("-")
        data["link"].append('<a href= http://crispr.i2bc.paris-saclay.fr/'
                            'crispr/crispr_db.php?checked%5B%5D=NC_007795 '
                            'target="_blank">CRIPSRdb</a>')
        check_tran(crispr, gffs["transcript"], data)
        check_tss(crispr, gffs["TSS"], "tss", data, "CRISPR")
        check_ps(crispr, gffs["processing_site"], "ps", data)
        check_term(crispr, gffs["terminator"], 30, 30, data)
        get_quant(crispr, quants, data, "CRISPR")

def get_srna_info(srnas, gffs, data, quants):
    for srna in srnas:
        import_basic(data, srna, "ncRNA")
        if "sRNA_0" not in srna.attributes["Name"]:
            data["gene_name"].append(srna.attributes["Name"])
        else:
            data["gene_name"].append("-")
        filename = "_".join([srna.feature, str(srna.start),
                             str(srna.end), srna.strand]) + ".html"
        data["link"].append('<a href=' + filename +
                            ' target="_blank">Co-epxression analysis</a>')
        check_tran(srna, gffs["transcript"], data)
        check_tss(srna, gffs["TSS"], "tss", data, "ncRNA")
        check_ps(srna, gffs["processing_site"], "ps", data)
        check_term(srna, gffs["terminator"], 30, 30, data)
        get_quant(srna, quants, data, "ncRNA")

def get_sorf_info(sorfs, gffs, data, quants):
    for sorf in sorfs:
        import_basic(data, sorf, "sORF")
        data["gene_name"].append("-")
        data["link"].append("NA")
        check_tran(sorf, gffs["transcript"], data)
        check_tss(sorf, gffs["TSS"], "tss", data, "sORF")
        check_ps(sorf, gffs["processing_site"], "ps", data)
        check_term(sorf, gffs["terminator"], 30, 300, data)
        get_quant(sorf, quants, data, "sORF")

def gen_column():
    columns = [
            TableColumn(
                field="features", title="Features", width=300,
                formatter=HTMLTemplateFormatter(
                template='<a href=<%= value %>.html target="_blank"><%= value %></a>')),
            TableColumn(field="start", title="Begin", width=150),
            TableColumn(field="end", title="End", width=150),
            TableColumn(field="strand", title="Strand", width=100),
            TableColumn(field="gene_name", title="Gene name", width=400),
            TableColumn(field="tran", title="Parent transcript"),
            TableColumn(field="tss", title="Associated TSSs", width=1000),
            TableColumn(field="ps", title="Associated processing sites", width=800),
            TableColumn(field="term", title="Associated terminators", width=420),
            TableColumn(
                field="link", title="Functional reference", width=460,
                formatter=HTMLTemplateFormatter(template='<%= value %>')),
            TableColumn(field="TSB_OD_0.2", title="RPKM:TSB_OD_0.2"),
            TableColumn(field="TSB_OD_0.2_TEX", title="RPKM:TSB_OD_0.2_TEX"),
            TableColumn(field="TSB_OD_0.5", title="RPKM:TSB_OD_0.5"),
            TableColumn(field="TSB_OD_0.5_TEX", title="RPKM:TSB_OD_0.5_TEX"),
            TableColumn(field="TSB_OD_1", title="RPKM:TSB_OD_1"),
            TableColumn(field="TSB_OD_1_TEX", title="RPKM:TSB_OD_1_TEX"),
            TableColumn(field="TSB_t0", title="RPKM:TSB_t0"),
            TableColumn(field="TSB_t0_TEX", title="RPKM:TSB_t0_TEX"),
            TableColumn(field="TSB_t1", title="RPKM:TSB_t1"),
            TableColumn(field="TSB_t1_TEX", title="RPKM:TSB_t1_TEX"),
            TableColumn(field="TSB_t2", title="RPKM:TSB_t2"),
            TableColumn(field="TSB_t2_TEX", title="RPKM:TSB_t2_TEX"),
            TableColumn(field="TSB_ON", title="RPKM:TSB_ON"),
            TableColumn(field="TSB_ON_TEX", title="RPKM:TSB_ON_TEX"),
            TableColumn(field="pMEM_OD_0.2", title="RPKM:pMEM_OD_0.2"),
            TableColumn(field="pMEM_OD_0.2_TEX", title="RPKM:pMEM_OD_0.2_TEX"),
            TableColumn(field="pMEM_OD_0.5", title="RPKM:pMEM_OD_0.5"),
            TableColumn(field="pMEM_OD_0.5_TEX", title="RPKM:pMEM_OD_0.5_TEX"),
            TableColumn(field="pMEM_OD_1", title="RPKM:pMEM_OD_1"),
            TableColumn(field="pMEM_OD_1_TEX", title="RPKM:pMEM_OD_1_TEX"),
            TableColumn(field="pMEM_t0", title="RPKM:pMEM_t0"),
            TableColumn(field="pMEM_t0_TEX", title="RPKM:pMEM_t0_TEX"),
            TableColumn(field="pMEM_t1", title="RPKM:pMEM_t1"),
            TableColumn(field="pMEM_t1_TEX", title="RPKM:pMEM_t1_TEX"),
            TableColumn(field="pMEM_t2", title="RPKM:pMEM_t2"),
            TableColumn(field="pMEM_t2_TEX", title="RPKM:pMEM_t2_TEX"),
            TableColumn(field="pMEM_ON", title="RPKM:pMEM_ON"),
            TableColumn(field="pMEM_ON_TEX", title="RPKM:pMEM_ON_TEX"),
            TableColumn(field="frag", title="Fragmentation")]
    return columns

def gen_html(filename):
    if args.features != "all":
        output_file(filename, mode='inline')

def main():
    gffs = {}
    names = {}
    quants = []
    if args.features == "all":
        output_file("main.html", mode='inline')
    data = {"link": [], "features": [], "start": [], "end": [], "strand": [],
            "gene_name": [], "pMEM_OD_0.2": [], "pMEM_OD_0.5": [],
            "pMEM_OD_1": [], "pMEM_t0": [], "pMEM_t1": [], "pMEM_t2": [],
            "pMEM_ON": [], "TSB_OD_0.2": [], "TSB_OD_0.5": [],
            "TSB_OD_1": [], "TSB_t0": [], "TSB_t1": [], "TSB_t2": [],
            "TSB_ON": [], "pMEM_OD_0.2_TEX": [], "pMEM_OD_0.5_TEX": [],
            "pMEM_OD_1_TEX": [], "pMEM_t0_TEX": [], "pMEM_t1_TEX": [],
            "pMEM_t2_TEX": [], "pMEM_ON_TEX": [], "TSB_OD_0.2_TEX": [],
            "TSB_OD_0.5_TEX": [], "TSB_OD_1_TEX": [], "TSB_t0_TEX": [],
            "TSB_t1_TEX": [], "TSB_t2_TEX": [], "TSB_ON_TEX": [],
            "frag": [], "tss": [], "ps": [], "term": [], "tran": []}
    nh = open(args.name_file, "r")
    for row in csv.reader(nh, delimiter='\t'):
        names[row[0]] = row[3]
    fh = open(args.gene_quanti_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if not row[0].startswith("Orientation"):
            quants.append(row)
    gff_f = open(args.input_file, "r")
    for entry in Gff3Parser().entries(gff_f):
        if entry.feature not in gffs.keys():
            gffs[entry.feature] = []
        gffs[entry.feature].append(entry)
    for feature in gffs.keys():
        if (feature == args.features) or (args.features == "all"):
            if feature == "riboswitch":
                gen_html("riboswitch.html")
                get_rfam_info(gffs["riboswitch"], gffs, data, quants,
                              "riboswitch")
            elif feature == "RNA_thermometer":
                gen_html("RNA_thermometer.html")
                get_rfam_info(gffs["RNA_thermometer"], gffs, data,
                              quants, "RNA_thermometer")
            elif feature == "CRISPR":
                gen_html("CRISPR.html")
                get_crispr_info(gffs["CRISPR"], gffs, data, quants)
            elif feature == "tRNA":
                gen_html("tRNA.html")
                get_rna_info(gffs["tRNA"], gffs, data, quants, "tRNA")
            elif feature == "rRNA":
                gen_html("rRNA.html")
                get_rna_info(gffs["rRNA"], gffs, data, quants, "rRNA")
            elif feature == "sORF":
                gen_html("sORF.html")
                get_sorf_info(gffs["sORF"], gffs, data, quants)
            elif feature == "ncRNA":
                gen_html("ncRNA.html")
                get_srna_info(gffs["ncRNA"], gffs, data, quants)
            elif feature == "CDS":
                gen_html("CDS.html")
                get_cds_info(gffs["CDS"], gffs, data, quants, names)
    source = ColumnDataSource(data)
    columns = gen_column()
    data_table = DataTable(source=source, columns=columns, height=3000, width=6000)
    save(widgetbox(data_table))
    
if __name__ == "__main__":
    main()
