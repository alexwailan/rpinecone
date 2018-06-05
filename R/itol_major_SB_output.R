#' Colourstrip Itol Output Major sub-lineages.
#'
#' Function to output the results of pinecone to itol
#' @param input output from pinecone
#' @import RColorBrewer
#' @export

itol_major_SL_output <- function(input){

  majSLno <- input$majorSLno
  tips <- input$ntips
  table <- input$table

  Date <- Sys.Date ()
  colnames (table) <- NULL


  n <- majSLno
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

  #Itol Legend vectors
  LEGEND_SHAPES_vector <- paste(rep(1,n), collapse=", ")
  LEGEND_COLORS_vector <- paste(colour_list[1:n], collapse=",")
  LEGEND_LABELS_vector <- paste(seq(1:n), collapse=",")

  #Itol Header Settings
  outputTableHeader <- rbind("DATASET_COLORSTRIP",
                             "SEPARATOR COMMA",
                             "DATASET_LABEL,Major Sub-lineage",
                             "COLOR_BRANCHES,0",
                             "LEGEND_TITLE,Major Sub-lineage",
                             paste("LEGEND_SHAPES,",LEGEND_SHAPES_vector, sep=""),
                             paste("LEGEND_COLORS,",LEGEND_COLORS_vector, sep=""),
                             paste("LEGEND_LABELS,",LEGEND_LABELS_vector, sep=""),
                             "DATA")
  for (i in 1:9){ #Place headers into Table
    outputTable[i] <- outputTableHeader[i, ]
  }

  # Collapse each row of the above tableinto one element for export after headers
  for (i in 1:nrow(preparingTable)){
    row = i + 9
    outputTable[row] <- paste(preparingTable[i,],collapse=",")
  }

  outputname <- paste("major_subleinages_itol_output_", Date, ".txt", sep = "")
  write.table(outputTable,file=outputname,sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

}
