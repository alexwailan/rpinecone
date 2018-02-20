#' Colourstrip Itol Output Major sub-groups.
#'
#' Function to output the results of maria to itol
#' This will change your tip labels to the sub-group number of the isolate
#' @param input output from maria
#' @import RColorBrewer

itol_major_SB_output <- function(input){

  majSBno <- input$majorSBno
  tips <- input$ntips
  table <- input$output

  Date <- Sys.Date ()
  colnames (table) <- NULL


  n <- majSBno
  preparingTable = mat.or.vec(tips, 3)

  for(i in 1:tips){
    preparingTable[i, ] <- table[i, c(1, 3, 3)] #Take the samples and their major cluster dupicated in two columns
  }

  for(j in 1:n){ #If less than 12 Major Groups - Generate Qualitative Colours and replace the column 2 with a colour for each Major Sub-group
    brew_list <- brewer.pal(12, "Paired")
    colour_list <- colorRampPalette(brewer.pal(12, "Paired"))(n)
    preparingTable[, 2][preparingTable[, 2]==j] <- colour_list[j] #For said number, check in the second column, which element said number is equal to; for those elements replace with colour HEX code
  }

  preparingTable <- preparingTable[-grep('0', preparingTable[, 2]), ]

  outputTable = mat.or.vec(nrow(preparingTable) + 9, 1)

  #Itol Header Settings
  outputTableHeader <- rbind("DATASET_COLORSTRIP",
                                   "SEPARATOR COMMA",
                                   "DATASET_LABEL,Major Sub-Group",
                                   "COLOR_BRANCHES,0",
                                   "LEGEND_TITLE,Major Sub-Group",
                                   "LEGEND_SHAPES",
                                   "LEGEND_COLORS",
                                   "LEGEND_LABELS",
                                   "DATA")
  for (i in 1:9){ #Place headers into Table
    outputTable[i] <- outputTableHeader[i, ]
  }

  # Collapse each row of the above tableinto one element for export after headers
  for (i in 1:nrow(preparingTable)){
    row = i + 9
    outputTable[row] <- paste(preparingTable[i,],collapse=",")
  }

  outputname <- paste("major_subgroups_itol_output_", Date, ".txt", sep = "")
  write.table(outputTable,file=outputname,sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

}
