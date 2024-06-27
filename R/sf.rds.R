sf.rds <- function(variable,file = "viarablename"){
  dir.create(file)
  namefile <- gsub(" ","",paste0("CWAS_result",file,Sys.time(),".rds"))
  namefile <- gsub(":","",namefile)
  #print(variable)
  saveRDS(variable,file = paste0(file,"/",namefile),compress = F)
  output <- paste0("CWAS_result ",namefile,"Has been saved in",file,"(Folder)")
  return(output)
}


