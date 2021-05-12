 args <- commandArgs(TRUE)
input <- args[1]
ouput <- gsub("indel|INDEL", "indel_hgvs", input)
in_df <- read.csv(input, header=TRUE, stringsAsFactors = FALSE)

 library(RCurl)
library(XML)
library(rjson)
 
 amino_df <- read.table("/haplox/tools/mhgvs/Amino_acid_descriptions_tab.txt",header=TRUE, stringsAsFactors = FALSE,sep="\t")
 
  amino_df_subset <- amino_df[-c(22,25),1:2]
 
  myGet <- function(x){
	tmp <- gsub("\\(", "", unlist(strsplit(x, ":"))[2])
	tmp <- gsub("\\)", "", tmp)
	return(tmp)
  }
 
 
for(j in seq(nrow(in_df))){
	print(j)
	print(in_df[j,])
if ( (!grepl("->", in_df[j, "base"])) && grepl("p\\.", in_df[j, "AA"])){
	if(in_df[j, "NM"] != "ARID1B"){
		
		input <- paste0("https://mutalyzer.nl/json/runMutalyzerLight?variant=", in_df[j, "NM"],":",in_df[j, "base"])
		
	}else if(in_df[j, "NM"] == "ARID1B"){
		input <- paste0("https://mutalyzer.nl/json/runMutalyzerLight?variant=", "NM_017519.2",":",in_df[j, "base"])
		print("ARID1B NM_017519.2")
	}
#    datas <- getURL("https://mutalyzer.nl/json/runMutalyzerLight?variant=NM_005228:c.2235_2249del")
	datas <- getURL(input)
	datas_json <-  fromJSON(datas)
#   	
	res_vec2 <- gsub("NM.*:(.*)", "\\1", datas_json$proteinDescriptions)  #[1] "NM_005228(EGFR_i001):p.(Glu746_Ala750del)"
	res_vec2 <- gsub("\\)", "", gsub("\\(", "", res_vec2))
	res_vec1 <- gsub("NM.*:(.*)", "\\1", datas_json$transcriptDescriptions) #[1] "NM_005228(EGFR_v001):c.2235_2249del"
  ### remove X, *
#  sapply(contents, myGet)
   for(i in seq(nrow(amino_df_subset))){
	res_vec2 <- gsub(amino_df_subset[i, 2], amino_df_subset[i, 1], res_vec2)
   }
   in_df[j, "base"] <-  res_vec1
   in_df[j, "AA"]   <-  res_vec2
}#if ( (!grepl("->", in_df[j, "base"])) && grepl("p\\.", in_df[j, "AA"])){

}  
write.table(in_df, file=ouput, row.names = FALSE, na = "", quote=FALSE, sep=",")

