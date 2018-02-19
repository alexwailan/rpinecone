#' itol_labels_template
#'
#' Function to output the results of maria to itol
#' This will change your tip labels to the sub-group number of the isolate
#' @param input output from maria
#'

itol_labels_template <- function(input){

  tips <- input$ntips
  table <- input$itolOutput
  Date <- Sys.Date()
  colnames (table) <- NULL
  outputTable = mat.or.vec(tips +3, 1)

  #Itol Header Settings
  outputTableHeaderVector <- rbind ("LABELS",
                                   "SEPARATOR COMMA",
                                   "DATA"
  )

  for (i in 1:3){ #Place headers into Table
    outputTable[i] <- outputTableHeaderVector[i,]
  }
  preparingTable = mat.or.vec(tips,2)

  for (i in 1:tips){

    #Take the samples and their major cluster dupicated in two columns
    preparingTable[i,] <- table[i,c(1,2)]

  }

  #With rows of the above table, collapse each row and its elements into one element for export after headers
  for (i in 1:nrow(table)){
    row = i +3
    outputTable[row] <- paste(preparingTable[i,],collapse=",")
  }

  subgroupoutput <- paste("subgroups_itol_label_output_", Date, ".txt", sep = "")
  write.table(outputTable, file = subgroupoutput, sep = ",", col.names = F, row.names = F, quote = F)

}
