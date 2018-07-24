#scMethylation_script.R

library(Biostrings)
library(bsseq)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(readr)
library(HDF5Array)

#  Replace  this  with  the  location/names  of  the  output  from  bismark_methylation_extractor
covgz_files  <-  c("../data/cell_1.bismark.cov.gz",
                   "../data/cell_2.bismark.cov.gz",                                  
                   "../data/cell_3.bismark.cov.gz")

getMethCov <- function(covgz_file, gr) {
  tab <- read_tsv(covgz_file,
                  col_types = "ciiiii",
                  col_names=c("chr", "pos", "pos2", "meth_percent", "m_count", "u_count"))
  tab_gr <- GRanges(tab$chr, IRanges(tab$pos, tab$pos))
  m <- rep(0, length(gr))
  cov <- rep(0, length(gr))
  ovl <- suppressWarnings(findOverlaps(tab_gr, gr))
  m[subjectHits(ovl)]   <- tab$m_count[queryHits(ovl)]
  cov[subjectHits(ovl)] <- tab$m_count[queryHits(ovl)] + tab$u_count[queryHits(ovl)]
  return(list(m=m, cov=as.integer(cov)))
}

#  Set  up  genome-wide  CpG  GRanges
#  On  the  plus  strand  we  keep  the  left-most  position  of  the  match
#  On  the  minus  strand  we  keep  theright-most  position  of  the  match

cpg_gr  <-  DNAString("CG")
cpg_gr  <-  vmatchPattern(cpg_gr,  BSgenome.Hsapiens.NCBI.GRCh38)
cpg_gr  <-  keepStandardChromosomes(cpg_gr, pruning.mode="coarse")

s  <-  start(cpg_gr)
e  <-  end(cpg_gr)
plus_idx  <-  as.logical(strand(cpg_gr)=="+")
minus_idx  <-  as.logical(strand(cpg_gr)=="-")
e[plus_idx]  <-  s[plus_idx]    #  Plus  strand
s[minus_idx]  <-  e[minus_idx]  #  Minus  strand
start(cpg_gr)  <-  s  
end(cpg_gr)  <-  e 

#  Get  methylation,  coverage  matrices  and  store  on  disk  as  HDF5Arrays
hdf5_m  <-  list()


hdf5_cov  <-  list()
for(covgz_file  in  covgz_files)  {    
  tmp  <-  getMethCov(covgz_file,  cpg_gr)    
  m  <-  tmp$m    
  cov  <-  tmp$cov    
  samplename  <-  basename(covgz_file)    
  samplename  <-  gsub(samplename, pattern="_PE_report.txt", replacement="")    
  hdf5_file  <-  paste0(samplename,  ".hdf5")   
  if(file.exists(hdf5_file)) file.remove(hdf5_file)    
  hdf5_m[[samplename]]  <-  writeHDF5Array(matrix(m),  name="m",  file=hdf5_file)    
  hdf5_cov[[samplename]]  <-  writeHDF5Array(matrix(cov),  name="cov",  file=hdf5_file)
}

M  <-  do.call("cbind", hdf5_m)
Cov  <-  do.call("cbind", hdf5_cov)

#  Create  a  bsseq  object  ready  for  downstream  analysis
bs  <-  BSseq(gr=cpg_gr,  M=M,  Cov=Cov)

