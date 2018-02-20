#' Colourstrip Itol Output of Sub-groups.
#'
#' Function to output the results of maria to itol.
#' This will change your tip labels to the sub-group number of the isolate.
#' @param input Output list from main maria function.
#' @import RColorBrewer

itol_subGroup_output <- function(input){

  tips <- input$ntips
  table <- input$output
  SGNo <- input$SGNo

  Date <- Sys.Date ()
  colnames (table) <- NULL

  preparingTable = mat.or.vec(tips, 3)

  for(i in 1:tips){
    preparingTable[i, ] <- table[i, c(1, 2, 2)] #Take the samples and their major cluster dupicated in two columns
  }

  for(l in 1:SGNo){
    brew_list <- brewer.pal(12, "Paired")
    colour_list <- colorRampPalette(brewer.pal(12, "Paired"))(SGNo)
    preparingTable[, 2][preparingTable[, 2]==l] <- colour_list[l] #For said number, check in the second column, which element said number is equal to; for those elements replace with colour HEX code
  }

  preparingTable <- preparingTable[-grep('singleton', preparingTable[, 2]), ]

  outputTable = mat.or.vec(nrow(preparingTable) + 9, 1)

  #Itol Header Settings
  outputTableHeaderVector <- rbind("DATASET_COLORSTRIP",
                                   "SEPARATOR COMMA",
                                   "DATASET_LABEL,Sub-Group",
                                   "COLOR_BRANCHES,0",
                                   "LEGEND_TITLE,Sub-Group",
                                   "LEGEND_SHAPES",
                                   "LEGEND_COLORS",
                                   "LEGEND_LABELS",
                                   "DATA")
  for (i in 1:9){ #Place headers into Table
    outputTable[i] <- outputTableHeaderVector[i, ]
  }


  # Collapse each row of the above tableinto one element for export after headers
  for (i in 1:nrow(preparingTable)){
    row = i + 9
    outputTable[row] <- paste(preparingTable[i,],collapse=",")
  }

  outputname <- paste("subgroups_itol_output_", Date, ".txt", sep = "")
  write.table(outputTable,file=outputname,sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

}
