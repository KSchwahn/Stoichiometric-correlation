
write_list <-function(my_list, wb_name = 'var1.xlsx') { 
#funtion to write R-list into Microsoft-Excel Workbook
  #each entry of the list will become a seperated worksheet
  #my_list: contains list of R objects (tested with vector, matrix and data.frame)
  #if names(my_list) is set, the workseet will be named
  #wb_name: name of the workbook created
  #the XLConnect library needs Java to be installed.
  
  #if no XLConnect package is installed .csv files are generated instead.
  #each list entry will be written into a separated file
  #the file names are automatically adjusted
  
  test_var =  "XLConnect" %in% rownames(installed.packages())
  
  if(test_var == T){
    library(XLConnect)
    wb <- loadWorkbook(wb_name, create = TRUE)
    createSheet(wb, names(my_list))
    writeWorksheet(wb, my_list, names(my_list),header=T)
    saveWorkbook(wb)
  }
  
 else {
    name = gsub("(.*)\\..*[a-z]", "\\1",wb_name)
    for(i in 1:length(my_list)){
      file_name = paste(name,"_",i,".csv",sep="")
      write.table(x=my_list[[i]],file=file_name,quote = F,col.names = T,row.names = F,sep=",")
    }
  }
}
