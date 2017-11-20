main(){
#    correlation_one
    correlation_all
#    network
#    gen_index
}

correlation_one(){
#    mkdir query
#    python3 bin/corr_query_old.py -pc 0.77 -nc -0.77 -p 1832755:1832968:+ -g input/GeneSpecificInformation_NCTC8325_del.csv -i output/gene_wise_quantifications_combined_deseq2_new.csv
     python3 bin/corr_query_one.py -go bin/goatools/data/go_database \
                               -gp bin/goatools/scripts/find_enrichment.py \
                               -ga bin/goatools/data/go_database \
                               -po bin/goatools/data/group_pop \
                               -gb input/go-basic.obo \
                               -pc 0.77 -nc -0.77 -p 2781669:2781891:- \
                               -n input/GeneSpecificInformation_NCTC8325_del.csv \
                               -g input/Staphylococcus_aureus_HG003_merge_features.gff \
                               -s Staphylococcus_aureus_HG003 \
                               -i input/gene_wise_quantifications_combined_deseq2_final_poly_u.csv > test1
#     mv *.html html
#    mv *_positive* query
#    mv *_negative* query
}

correlation_all(){
#    mkdir query
    python3 bin/corr_query.py -go bin/goatools/data/go_database \
                              -gp bin/goatools/scripts/find_enrichment.py \
                              -ga bin/goatools/data/go_database \
                              -po bin/goatools/data/group_pop \
                              -gb input/go-basic.obo \
                              -pc 0.77 -nc -0.77 \
                              -g input/Staphylococcus_aureus_HG003_merge_features.gff \
                              -n input/GeneSpecificInformation_NCTC8325_del.csv \
                              -s Staphylococcus_aureus_HG003 \
                              -i input/gene_wise_quantifications_combined_deseq2_final_poly_u.csv
    mv *.html html
}

network(){
    python3 bin/network_bokeh.py -i output/gene_wise_quantifications_combined_deseq2_final_poly_u.csv
#    python3 bin/network_bokeh2.py -i output/gene_wise_quantifications_combined_deseq2_final_poly_u.csv
}

gen_index(){
    for FEATURE in riboswitch RNA_thermometer CDS tRNA rRNA sORF ncRNA CRISPR all
    do 
        python3 bin/table.py -i input/Staphylococcus_aureus_HG003_merge_features.gff \
                             -g input/gene_wise_quantifications_combined_rpkm.csv \
                             -n input/GeneSpecificInformation_NCTC8325_del.csv \
                             -f ${FEATURE}
    done
    mv *.html html
}

main
