#function section special for these investigation
write_list <-function(my_list, wb_name = 'var1.xlsx') { 
  library(XLConnect)
  wb <- loadWorkbook(wb_name, create = TRUE)
  createSheet(wb, names(my_list))
  writeWorksheet(wb, my_list, names(my_list),header=T)
  saveWorkbook(wb)
}
