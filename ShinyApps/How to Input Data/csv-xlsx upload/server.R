if(!require(readxl)) install.packages(readxl)
if(!require(shiny)) install.packages(shiny)
if(!require(ggpubr)) install.packages(ggpubr)
if(!require(dplyr)) install.packages(dplyr)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 90*1024^2)

function(input, output) {
  Dataset <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    data = read_excel(path = inFile$datapath,col_names = FALSE, sheet = "Scenarios" )
    dims =  dim(data)
    pointer = as.data.frame(data[1:3,])
    
    pointer = as.data.frame(apply(pointer,2,function(x)  round((as.numeric(x)),2)))
    
    which.is = c(input$PLND_sens,
                 input$New_Sens,
                 input$New_Spes)
    
    col = which(apply(pointer, 2, function(x) identical(x,which.is)))
    
    
    data = data[(7:dims[1]), (col:(col+3))]
    data = as.data.frame(apply(data,2, as.numeric ))
    names(data)= c("Nano_incremental_costs","Nano_incremental_QALYs","PSMA_PET_CT_incremental_costs","PSMA_PET_CT_incremental_QALYs")
    return(data)
  })
  
  # Show table:
  output$table <- renderTable({
    
    Dataset2 = head(Dataset(), n = 10)
    
    return(Dataset2)
  })
  
  # Show plot1:
  output$plot1 <- renderPlot({
    data = Dataset()
    plot1 = ggscatter(data = data,"Nano_incremental_QALYs", "Nano_incremental_costs") + 
      geom_point(data = data[1,],colour = "red" )+  # this adds a red point
      geom_text(data=data[1,], label="point estimates nano-MRI", vjust=1,hjust=1.1) # this adds a label for the red point

    
    return(plot1)
  })
  
  # Show plot2:
  output$plot2 <- renderPlot({
    data = Dataset()
    plot2 = ggscatter(data = data,"PSMA_PET_CT_incremental_QALYs", "PSMA_PET_CT_incremental_costs") + 
      geom_point(data = data[1,],colour = "red" )+  # this adds a red point
      geom_text(data=data[1,], label="point estimates nano-MRI", vjust=1,hjust=1.1) # this adds a label for the red point
    
    return(plot2)
  })
  
}