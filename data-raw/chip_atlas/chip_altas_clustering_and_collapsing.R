dgCMatrix_to_01_data_frame = function(dgCM, cutoff = 100, trim = FALSE) {
  dgCM[ dgCM <= cutoff ] = 0
  dgCM[ dgCM >  cutoff ] = 1
  dgCM = data.frame( data.matrix(dgCM) )
  if(trim) {
    colnames(dgCM) = sub( "_SRX.*" ,"", colnames(dgCM))
    colnames(dgCM) = sub( "_Pluripotent.stem.cell_hESC" ,"_hESC", colnames(dgCM))
  }
  dgCM
}
collapse_tf_signatures = function( input_mat, input_annot){
  tf_mat = matrix(0, nrow = nrow( input_mat ), ncol = length(unique(input_annot$TF)))
  colnames(tf_mat) = unique(input_annot$TF)
  rownames(tf_mat) = rownames(input_mat)
  for ( j in 1: ncol (tf_mat)){
    mat_sel_tf = chip_atlas_symbol[  ,grep( colnames(tf_mat)[j] , colnames(input_mat))]
    mat_sel_tf = dgCMatrix_to_01_data_frame ( mat_sel_tf ,trim = TRUE)
    tf_peak_count_per_gene = rowSums(mat_sel_tf)
    FT_sig  = tf_peak_count_per_gene > max(tf_peak_count_per_gene)/4
    tf_mat[FT_sig ,j] = 1
  }
  tf_mat
}


data("chip_atlas_symbol", package = "xcore")
data("chip_atlas_meta", package = "xcore")

chip_atlas_meta = chip_atlas_meta[ !( chip_atlas_meta$TF %in% c("Epitope tags", "GFP", "Hepatitis B Virus X antigen" , "Biotin","Cyclobutane pyrimidine dimers",
                                                                "H3K9K14K18K23K27ac", "Pan-acetyllysine", "8-Hydroxydeoxyguanosine") )  ,]



tb = table( chip_atlas_meta$TF )
singos = chip_atlas_meta [ chip_atlas_meta$TF %in% names(tb)[tb == 1 ] , ]
duos   = chip_atlas_meta [ chip_atlas_meta$TF %in% names(tb)[tb == 2 ] , ]
fewos  = chip_atlas_meta [ chip_atlas_meta$TF %in% names(tb)[tb > 2 & tb < 10 ] , ]

tf_X_origin = table( muchos$TF , muchos$origin)
tf_X_origin [ (rowSums(tf_X_origin > 2) ) > 5 , ]


muchos = chip_atlas_meta [ chip_atlas_meta$TF %in% names(tb)[tb >= 10] , ]
muchos_mat = chip_atlas_symbol [ , colnames(chip_atlas_symbol) %in% muchos$name ]
table( colnames(muchos_mat) == muchos$name)

muchos_mat = collapse_tf_signatures ( input_mat = muchos_mat, input_annot = muchos)

muchos_mat = muchos_mat[ , colSums(muchos_mat) > 30]

save(muchos_mat, file = "data-raw/chip_atlas/muchos_mat.rda")

tail( sort()) )
\
tail( sort( table ( chip_atlas_meta$TF)))

tb = table( paste( chip_atlas_meta$TF , chip_atlas_meta$origin))


tb = table( chip_atlas_meta$TF)
tb = tb [ tb > 9]



chip_atlas_meta_sel = chip_atlas_meta [ chip_atlas_meta$TF %in% names(tb) ,]




tf_mat = as.data.frame(tf_mat)
rownames(tf_mat)[ tf_mat$TP53 == 1]
sort(colSums(tf_mat))

hist(colSums(tf_mat), breaks = 100)


hist(MYC_count , breaks = 30)


mat = MYC

FT_signature           = chip_peak_count_per_gene > max(chip_peak_count_per_gene)/2

barplot(table(chip_peak_count_per_gene)[1:10])




mean(sort( colSums(MYC), decreasing = T))
[1] 3539.966
> sd(sort( colSums(MYC), decreasing = T))
[1] 4196.798

barplot( sort( colSums(MYC), decreasing = T), las = 2,
         ylab = "# of genes with ChipSeq peak", main = "MYC")

barplot( sort( rowSums(MYC), decreasing = T), las = 2,
         ylab = "# of genes ChipSeq peak", main = "MYC")


chip_atlas_meta$primary_tissue = gsub( "Primary Tissue=|Tissue=|\\|.*" , "", chip_atlas_meta$description )
chip_atlas_meta$diagnosis      = gsub( ".*Tissue Diagnosis=" , "", chip_atlas_meta$description )

chip_atlas_meta_normal = chip_atlas_meta [ grep( "Normal" , chip_atlas_meta$description),]

chip_atlas_meta_normal = chip_atlas_meta [ grep( "Normal" , chip_atlas_meta$description),]



chip_atlas_meta$primary_tissue = gsub( "Primary Tissue=|Tissue=|\\|.*" , "", chip_atlas_meta$description )
tail(sort(table(sxr_meta$V4)))

sxr_meta = read.delim("~/projects/resources/chip_atlas/experimentList.tab_TF_hg19.txt", sep ="\t", header = F, quote = "\"")



sxr_meta [ sxr_meta$V4 %in% c("Epitope tags", "GFP"),][1:12,]

dim( chip_atlas_meta [ grep( "MeSH",chip_atlas_meta$description ),] )

sxr_meta[ grep( "MeSH",sxr_meta$V7 ),][1:12,]
