
msigdb = xcore::read_msigdb_gmts(path_to_gmts = "/Users/bogumil/projects/resources/msigdb/v7.0/" ,
                                 msigdb_version = ".v7.0.entrez.gmt" )
table(sub( "__.*" ,"", names(msigdb) ))


msigdb_mat_sparse = Matrix::drop0( xcore::list_to_matrix(msigdb) )

msigdb_mat = msigdb_mat_sparse

usethis::use_data(msigdb_mat)
