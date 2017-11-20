#!/usr/bin/python

import os
import sys
import csv
import shutil
import argparse
import math
from subprocess import call
from scipy.stats.stats import pearsonr, spearmanr
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.spatial import distance
from bokeh.resources import INLINE 
from bokeh.io import output_file, save, hplot
from bokeh.layouts import widgetbox, column, row, layout, gridplot
from bokeh.models import Line, Label, FuncTickFormatter, ColumnDataSource, DataRange1d, Plot, LinearAxis, Grid, Circle, HoverTool, BoxSelectTool
from bokeh.models.widgets import (DataTable, TableColumn, DateFormatter,
                                  StringFormatter, NumberFormatter, HTMLTemplateFormatter,
                                  StringEditor, IntEditor, NumberEditor, SelectEditor)
from bokeh.sampledata.periodic_table import elements
from bokeh.plotting import figure
from bokeh.models.tickers import FixedTicker
from gff3 import Gff3Parser


__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--gene_quanti_file",help="gene quantification file")
parser.add_argument("-g","--gff_file", default=None, help="gene gff file")
parser.add_argument("-n","--name_file", default=None, help="gene name file")
parser.add_argument("-pc","--pos_cut", type=float, default=0.77, help="the cutoff of positive correlation, default is 0.8")
parser.add_argument("-nc","--neg_cut", type=float, default=-0.77, help="the cutoff of negative correlation, default is -0.8")
parser.add_argument("-go","--go_file", help="association go file")
parser.add_argument("-gp","--goatools_path", help="path of find_enrichment.py in goatools")
parser.add_argument("-gb","--obo_file", help="path fo go.obo")
parser.add_argument("-po","--population_file", help="Go population file")
parser.add_argument("-ga","--go_association", help="Go association file")
parser.add_argument("-s","--strain", help="Strain name")
parser.add_argument("-k","--known_srna_only", action="store_true", help="only analyze known sRNA")

args = parser.parse_args()

def get_quant(exp, members):
     members["TSB_OD_0.2"].append(float(exp[0]))
     members["TSB_OD_0.5"].append(float(exp[1]))
     members["TSB_OD_1"].append(float(exp[2]))
     members["TSB_ON"].append(float(exp[3]))
     members["TSB_t0"].append(float(exp[4]))
     members["TSB_t1"].append(float(exp[5]))
     members["TSB_t2"].append(float(exp[6]))
     members["pMEM_OD_0.2"].append(float(exp[7]))
     members["pMEM_OD_0.5"].append(float(exp[8]))
     members["pMEM_OD_1"].append(float(exp[9]))
     members["pMEM_ON"].append(float(exp[10]))
     members["pMEM_t0"].append(float(exp[11]))
     members["pMEM_t1"].append(float(exp[12]))
     members["pMEM_t2"].append(float(exp[13]))

def get_filename(info):
    attrs = "_".join(info.split("_")[4:])
    datas = attrs.split(";")
    for data in datas:
        if data.startswith("Name"):
            product = data.split("=")[-1]
            return product

def get_pro_id(info):
    attrs = "_".join(info.split("_")[3:])
    datas = attrs.split(";")
    for data in datas:
        if data.startswith("protein_id"):
            product = data.split("=")[-1]
            return product
    return None

def get_gene_name(genes, info, names):
    attrs = "_".join(info.split("_")[4:])
    datas = attrs.split(";")
    if info.split("_")[0] != "ncRNA":
        for data in datas:
            if data.startswith("locus_tag="):
                loc = data.split("=")[-1]
                break
        for locus, name in names.items():
            if locus == loc:
                if name != "-":
                    gene_name = name
                else:
                    gene_name = loc
        for data in datas:
            if data.startswith("Parent="):
                parent = data.split("=")[-1]
                break
        for gene in genes:
            if gene.attributes["ID"] in parent.split(","):
                link = "https://www.ncbi.nlm.nih.gov/gene/" + gene.attributes["db_xref"].split(":")[-1]
                hyper = '<a href= ' + link + ' target="_blank"> ' + gene.attributes["db_xref"] + '</a>'
        return hyper, gene_name
    else:
        gene_name = "novel"
        for data in datas:
            if data.startswith("Name") and ("sRNA_0" not in data):
                gene_name = data.split("=")[-1]
        filename = "_".join(info.split("_")[0:4]) + ".html"
        hyper = '<a href=' + filename + ' target="_blank">Co-epxression analysis</a>'
        return hyper, gene_name

def gen_input_goatools(filename, infos, exps, cutoff, q_info, q_exp, gos):
    out = open(filename, "w")
    datas = []
    ids = []
    pro_gos = []
    ccs = []
    for info in infos:
        if info != q_info:
            index = "_".join(info.split("_")[:4])
            index_q = "_".join(q_info.split("_")[:4])
            if "positive" in filename:
                if (float(spearmanr(exps[index], q_exp[index_q])[0]) >= cutoff):
                    if get_pro_id(info) is not None:
                        out.write(get_pro_id(info) + "\n")
                    detect = False
                    for pro_id, go_list in gos.items():
                        if pro_id in info:
                            pro_gos.append(go_list)
                            detect = True
                            break
                    if not detect:
                        pro_gos.append("NA")
                    datas.append(exps[index])
                    ids.append(info)
                    ccs.append("{0:.5f}".format(float(spearmanr(exps[index], q_exp[index_q])[0])))
            elif "negative" in filename:
                if (float(spearmanr(exps[index], q_exp[index_q])[0]) <= cutoff):
                    detect = False
                    if get_pro_id(info) is not None:
                        out.write(get_pro_id(info) + "\n")
                    for pro_id, go_list in gos.items():
                        if pro_id in info:
                            pro_gos.append(go_list)
                            detect = True
                            break
                    if not detect:
                        pro_gos.append("NA")
                    datas.append(exps[index])
                    ids.append(info)
                    ccs.append("{0:.5f}".format(float(spearmanr(exps[index], q_exp[index_q])[0])))
    out.close()
    out_go = open(filename + "_go", "w")
    call(["python3", args.goatools_path, "--pval=0.05", "--indent",
          "--obo=" + args.obo_file, filename, args.population_file,
          args.go_association], stdout=out_go)
    out_go.close()
    fh = open(filename + "_go", "r")
    start = False
    enrichs = []
    for row in csv.reader(fh, delimiter='\t'):
        if start:
            if row[2] == "e":
                enrichs.append(row[0].replace(".", ""))
        if row[0] == "GO":
            start = True
    fh.close()
    os.remove(filename)
    os.remove(filename + "_go")
    return datas, ids, pro_gos, enrichs, ccs

def import_basic(items, type_, go, cc, gene_name, link, exps, members):
    members["role"].append(type_)
    members["start"].append(items[1])
    members["end"].append(items[2])
    members["strand"].append(items[3])
    members["feature"].append(items[0])
    members["link"].append(link)
    members["gene_name"].append(gene_name)
    members["GO"].append(go)
    members["cc"].append(cc)
    for info, exp in exps.items():
        if info == "_".join(items[:4]):
            get_quant(exp, members)

def plot(filename, q_info, q_exp, names, pngname, infos, exps, cutoff, genes, gos, members, fig_title):
    tags = ["TSB_OD_0.2", "TSB_OD_0.5", "TSB_OD_1", "TSB_t0", "TSB_t1", "TSB_t2", "TSB_ON",
            "pMEM_OD_0.2", "pMEM_OD_0.5", "pMEM_OD_1", "pMEM_t0", "pMEM_t1", "pMEM_t2", "pMEM_ON"]
    datas, ids, pro_gos, enrichs, ccs = gen_input_goatools(
         filename, infos, exps, cutoff, q_info, q_exp, gos)
    f_datas = []
    f_ids = []
    f_ccs = []
    q_items = q_info.split("_")
    link, gene_name = get_gene_name(genes, q_info, names)
    if "positive" in filename:
        import_basic(q_items, "Query sRNA", "-", "-", gene_name, link, q_exp, members)
    for info, exp, go_list, cc in zip(ids, datas, pro_gos, ccs):
        link, gene_name = get_gene_name(genes, info, names)
        items = info.split("_")
        if go_list != "NA":
            for go in go_list:
                if go in enrichs:
                    f_datas.append(exp)
                    f_ids.append(info)
                    f_ccs.append(cc)
                    if "positive" in filename:
                        import_basic(items, "correlated", "\n".join(go_list), cc,
                                     gene_name, link, exps, members)
                    else:
                        import_basic(items, "anti-correlated", "\n".join(go_list), cc,
                                     gene_name, link, exps, members)
                    break
        else:
            f_datas.append(exp)
            f_ids.append(info)
            f_ccs.append(cc)
            if "positive" in filename:
                import_basic(items, "correlated", "-", cc,
                             gene_name, link, exps, members)
            else:
                import_basic(items, "anti-correlated", "-", cc,
                             gene_name, link, exps, members)
    hover = HoverTool(
        tooltips=[
            ("feature", "@feature"),
            ("start", "@start"),
            ("end", "@end"),
            ("strand", "@strand"),
            ("expression value", "@exp"),
            ("cc", "@cc")
        ]
    )
    x = np.arange(14)
    p = figure(title=fig_title, plot_width=850, plot_height=850,
               tools=["pan,wheel_zoom,box_zoom,reset,tap", hover]) ## tap for select line and highlight it
    for exp, id_, cc in zip(f_datas, f_ids, f_ccs):
        sources = {"feature": [], "start": [], "end": [], "strand": [], "exp": [], "cc": []}
        for e_v in exp:
            infos = id_.split("_")
            sources["feature"].append(infos[0])
            sources["start"].append(infos[1])
            sources["end"].append(infos[2])
            sources["strand"].append(infos[3])
            sources["exp"].append(e_v)
            sources["cc"].append(cc)
        l = p.line(x, exp, line_color="lightgrey", source=sources)
        l.selection_glyph = Line(line_color="blue", line_width=3)
        l.nonselection_glyph = Line(line_color="blue", line_width=2)
    for info, exp in q_exp.items():
        sources = {"feature": [], "start": [], "end": [], "strand": [], "exp": [], "cc": []}
        for e_v in exp:
            infos = info.split("_")
            sources["feature"].append("ncRNA")
            sources["start"].append(infos[1])
            sources["end"].append(infos[2])
            sources["strand"].append(infos[3])
            sources["exp"].append(e_v)
            sources["cc"].append("Query sRNA")
        l = p.line(x, exp, legend="Query sRNA", line_color="red", line_width=2, source=sources)
        l.selection_glyph = Line(line_color="red", line_width=3)
        l.nonselection_glyph = Line(line_color="red", line_width=2)
    p.xaxis.axis_label = "Conditions"
    p.yaxis.axis_label = "Log2 fold changes"
    p.xaxis.ticker = FixedTicker(ticks=x) ## set the number of ticks
    label_dict = {}
    for i, s in enumerate(tags):
        label_dict[i] = s
    p.xaxis.formatter = FuncTickFormatter(code="""
        var labels = %s;
        return labels[tick];
    """ % label_dict)  ## change the tick from number to string
    p.xaxis.major_label_orientation = math.pi/4 ## rotate labels of ticks
    label_opts = dict(
        x=670, y=0, text_font_size="8pt",
        x_units='screen', y_units='screen')
    msg1 = 'Highlighting the line by tapping it.'
    caption1 = Label(text=msg1, **label_opts)
    p.add_layout(caption1, 'above')
    return p

def gen_table(members):
    source = ColumnDataSource(members)
    columns = [
            TableColumn(field="role", title="Role", width=250),
            TableColumn(field="feature", title="Feature", width=250),
            TableColumn(field="start", title="Begin", width=120),
            TableColumn(field="end", title="End", width=120),
            TableColumn(field="strand", title="Strand", width=80),
            TableColumn(field="gene_name", title="Gene name", width=140),
            TableColumn(field="link", title="Reference link", formatter=HTMLTemplateFormatter(template='<%= value %>')),
            TableColumn(field="cc", title="C.C"),
            TableColumn(field="GO", title="GO term"),
            TableColumn(field="TSB_OD_0.2", title="Log2 fold change:TSB_OD_0.2"),
            TableColumn(field="TSB_OD_0.5", title="Log2 fold_change:TSB_OD_0.5"),
            TableColumn(field="TSB_OD_1", title="Log2 fold_change:TSB_OD_1"),
            TableColumn(field="TSB_t0", title="Log2 fold_change:TSB_t0"),
            TableColumn(field="TSB_t1", title="Log2 fold_change:TSB_t1"),
            TableColumn(field="TSB_t2", title="Log2 fold_change:TSB_t2"),
            TableColumn(field="TSB_ON", title="Log2 fold_change:TSB_ON"),
            TableColumn(field="pMEM_OD_0.2", title="Log2fold_change:pMEM_OD_0.2"),
            TableColumn(field="pMEM_OD_0.5", title="Log2 fold_change:pMEM_OD_0.5"),
            TableColumn(field="pMEM_OD_1", title="Log2 fold_change:pMEM_OD_1"),
            TableColumn(field="pMEM_t0", title="Log2 fold_change:pMEM_t0"),
            TableColumn(field="pMEM_t1", title="Log2 fold_change:pMEM_t1"),
            TableColumn(field="pMEM_t2", title="Log2 fold_change:pMEM_t2"),
            TableColumn(field="pMEM_ON", title="Log2 fold_change:pMEM_ON"),
        ]
    data_table = DataTable(source=source, columns=columns, editable=True, height=4000, width=5000)
    return data_table

def main():
    gos = {}
    goh = open(args.go_file, "r")
    for row in csv.reader(goh, delimiter='\t'):
        gos[row[0]] = row[1].split(";")
    genes = []
    names = {}
    srnas = []
    nh = open(args.name_file, "r")
    for row in csv.reader(nh, delimiter='\t'):
        names[row[0]] = row[3]
    gh = open(args.gff_file, "r")
    for entry in Gff3Parser().entries(gh):
        if entry.feature == "gene":
            genes.append(entry)
        elif entry.feature == "ncRNA":
            srnas.append(entry)
    for srna in srnas:
        members = {"role": [], "feature": [], "start": [], "end": [], "strand": [], "link": [],
                   "gene_name": [], "GO": [], "cc": [], "pMEM_OD_0.2": [], "pMEM_OD_0.5": [],
                   "pMEM_OD_1": [], "pMEM_t0": [], "pMEM_t1": [], "pMEM_t2": [],
                   "pMEM_ON": [], "TSB_OD_0.2": [], "TSB_OD_0.5": [],
                   "TSB_OD_1": [], "TSB_t0": [], "TSB_t1": [], "TSB_t2": [],
                   "TSB_ON": []}
        exps ={}
        infos = []
        print("_".join([srna.feature, str(srna.start), str(srna.end), srna.strand]))
        output_file("_".join([srna.feature, str(srna.start), str(srna.end), srna.strand]) + ".html", mode='inline')
        fh = open(args.gene_quanti_file, "r")
        for row in csv.reader(fh, delimiter='\t'):
            if (row[0] != "Orientation") and (row[0] == "sense"):
                if (not args.known_srna_only) or (
                       (args.known_srna_only) and ("Name=sRNA" not in row[9])):
                    exps["_".join([row[3], row[4], row[5], row[7]])] = [
                                 float(row[10]), float(row[11]), float(row[12]),
                                 float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                                 float(row[17]), float(row[18]), float(row[19]),
                                 float(row[20]), float(row[21]), float(row[22]), float(row[23])]
                infos.append("_".join([row[3], row[4], row[5], row[7], row[9]]))
                if (row[4] == str(srna.start)) and (row[5] == str(srna.end)) and (row[7] == str(srna.strand)):
                    q_exp = {"_".join([row[3], row[4], row[5], row[7]]): [
                             float(row[10]), float(row[11]), float(row[12]),
                             float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                             float(row[17]), float(row[18]), float(row[19]),
                             float(row[20]), float(row[21]), float(row[22]), float(row[23])]}
                    q_info = "_".join([row[3], row[4], row[5], row[7], row[9]])
        name = get_filename(q_info)
        if "/" in name:
            name = "_".join([name.replace("/", "_or_"), str(srna.start), str(srna.end), str(srna.strand)])
        else:
            name = "_".join([name, str(srna.start), str(srna.end), str(srna.strand)])
        p1 = plot(name + "_positive", q_info, q_exp, names,
                 name + "_positive.png", infos, exps, args.pos_cut, genes, gos, members, "Correlation")
        p2 = plot(name + "_negative", q_info, q_exp, names,
                 name + "_negative.png", infos, exps, args.neg_cut, genes, gos, members, "Anti-correlation")
        grid = gridplot([p1, p2], ncols=2)
        table = gen_table(members)
        l = layout([[grid], [widgetbox(table)]])
        save(l)
        fh.close()

if __name__ == "__main__":
    main()
