# Title     : TODO
# Objective : TODO
# Created by: chema
# Created on: 2023-03-03

library(Rsubread)

flattenGTF <- function(GTFfile,GTF.featureType="exon",GTF.attrType="gene_id",method="merge"){
  .check_string_param(GTF.featureType,"GTF.featureType")
  .check_string_param(GTF.attrType,"GTF.attrType")
  .check_string_param(method,"method")

  method <- match.arg(method,c("merge","chop"))
  GTFfile <- .check_and_NormPath(GTFfile, mustWork=TRUE, opt="GTFfile")

  fout <- file.path(".",paste0(".Rsubread_flattenGTF_pid",Sys.getpid()))
  cmd <- paste("RflattenGTF","-a",GTFfile,sep=.R_param_splitor)
  cmd <- paste(cmd,"-g",GTF.attrType,sep=.R_param_splitor)
  cmd <- paste(cmd,"-t",GTF.featureType,sep=.R_param_splitor)
  cmd <- paste(cmd,"-o",fout,sep=.R_param_splitor)
  if(method == 'chop') cmd <- paste(cmd,"-C",sep=.R_param_splitor)

  n <- length(unlist(strsplit(cmd,.R_param_splitor)))
  C_args <- .C("R_flattenGTF_wrapper",n,cmd,PACKAGE="Rsubread")
  if(file.exists(fout)){
	z <- read.delim(fout,stringsAsFactors=FALSE,colClasses=c("character","character","integer","integer","character"))
	file.remove(fout)
	z
  }else{
	warning("No output was generated.")
	data.frame()
  }
}

flattenGTF(

    # basic input/output options
    "/shared/sequences/ixodes_ricinus/ncRNA/Consensus/Consensus_ncRNA_gffcompare_with_OGS0.1_s90.CDS.mergeLoci.gtf",
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",

    # the option specifying the merging algorithm
    method = "merge")

