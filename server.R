library(bio3d)
library(shiny)
library(DT)


shinyServer(function(input, output) {
  
  sliderValues <- reactive({
    
    data.frame(
      Name =  "Threshold",
      Value = as.character(input$decimal), 
      stringsAsFactors=FALSE)
  }) 
  
  output$values <- renderTable({
    sliderValues()
  })
  
  slider <- reactive({
    
      value <- as.numeric(as.character(input$decimal))
      return(value)
  }) 
  
  get_kinesin_id <- reactive({
    
    kinesin_ids <- input$text
    kinesin_ids <- strsplit(kinesin_ids, ",")
    return(kinesin_ids)
  })
  
  
  get_output_table <- function(){
    
    fasta = "D:/final_alignment-h4-5_tempus2_barry_version.fa"
    file_fasta <- read.fasta(fasta, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
    
    conservation <- conserv(x=file_fasta$ali, method="similarity", sub.matrix="bio3d")
    conserv_mat <- as.data.frame(conservation)
    conserv_mat <- rbind(c(""), conserv_mat)
    
    
    consensus_table <- read.table("D:/alignment_guido_ntFull.xls", header = FALSE, sep = "\t")
    residue_table <- read.table("D:/sse_residue_nos_guido.xls", header = FALSE, sep = "\t")
    conservation_table <- read.table("D:/sse_conserv_tab_guido.xls", header = FALSE, sep = "\t")

    datalist = list()
    input = get_kinesin_id()
    input <- unlist(input)
    len <- length(input)
    output_table <- NULL
    output_1 <- as.data.frame(t(consensus_table[2,]))
    rows <- nrow(output_1)
    rownames(output_1) <- seq(rows)
    colnames(output_1) <- "Regions"
    output_5 <- as.data.frame(t(conservation_table[2,]), stringsAsFactors = F)
    rownames(output_5) <- seq(rows)
    colnames(output_5) <- "Conservation_Score"
    value <- slider()
    for(h in 1:len){
      input_1 <- input[h]
      input_1 <- toupper(input_1)
      for(i in 2:113){
        if(input_1 == as.character(consensus_table[i,1])){
          output_2 <- as.data.frame(consensus_table[i,])
          output_2 <- t(output_2)
          output_2[output_2 == "-" ] <- ""
          output_3 <- as.data.frame(residue_table[i-1,])
          output_3 <- t(output_3)
          output_3[is.na(output_3)] <- ""
        }
      }
      for(i in 2:609){
        if(as.numeric(as.character(conservation_table[2,i])) > value){
          if(output_3[i]!= ""){
            output_3[i] <- paste(output_3[i],'*',sep="")
          }
        }
      }
      output_4 <- as.data.frame(paste(t(output_2),t(output_3),sep="."))
      colnames(output_4) <- paste("Residues",toupper(input[h]),sep=".")
      datalist[[h]] <- output_4

    }
    big_data = do.call(cbind, datalist)
    options(stringsAsFactors = FALSE)
    output_table <- cbind(output_1, big_data, output_5)
    output_table <- output_table[-1,]
    rownames(output_table) <- seq(nrow(output_table))

    #datatable to create color bars according to the conservation score    

#     DT::datatable(output_table) %>% 
#       formatStyle(
#         'Conservation_Score',
#         background = styleColorBar(as.numeric(as.character(output_table$Conservation_Score)), 'steelblue'),
#         backgroundSize = '100% 90%',
#         backgroundRepeat = 'no-repeat',
#         backgroundPosition = 'center'
#       )
    
    return(output_table)
  }

  output$output_table <- renderDataTable({
    
    output_table <- NULL
    output_table <- get_output_table()
    return(output_table)
  })


})

