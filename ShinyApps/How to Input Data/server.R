library(foreign)
library(readr)
library(haven)
library(readxl)
library(xlsx)

shinyServer(function(input, output) {
  
  #Importing Data Code#
  ######### 
  
  ### Argument names:
  ArgNames <- reactive({
    Names <- names(formals(input$readFunction)[-1])
    Names <- Names[Names!="..."]
    return(Names)
  })
  
  ### write Argument names:
  write.ArgNames <- reactive({
    write.Names <- names(formals(input$write.Function)[-1])
    write.Names <- write.Names[write.Names!="..."]
    return(write.Names)
  })
  
  # Argument selector:
  output$ArgSelect <- renderUI({
    if (length(ArgNames())==0) return(NULL)
    
    selectInput("arg","Argument:",ArgNames())
  })
  
  # Argument selector:
  output$write.ArgSelect <- renderUI({
    if (length(write.ArgNames())==0) return(NULL)
    
    selectInput("write.arg","Argument:",write.ArgNames())
  })
  
  ## Arg text field:
  output$ArgText <- renderUI({
    fun__arg <- paste0(input$readFunction,"__",input$arg)
    
    if (is.null(input$arg)) return(NULL)
    
    Defaults <- formals(input$readFunction)
    
    if (is.null(input[[fun__arg]]))
    {
      textInput(fun__arg, label = "Enter value:", value = deparse(Defaults[[input$arg]])) 
    } else {
      textInput(fun__arg, label = "Enter value:", value = input[[fun__arg]]) 
    }
  })
  
  ## write.Arg text field:
  output$write.ArgText <- renderUI({
    write.fun__arg <- paste0(input$write.Function,"__",input$write.arg)
    
    if (is.null(input$write.arg)) return(NULL)
    
    write.Defaults <- formals(input$write.Function)
    
    if (is.null(input[[write.fun__arg]]))
    {
      textInput(write.fun__arg, label = "Enter value:", value = deparse(write.Defaults[[input$write.arg]])) 
    } else {
      textInput(write.fun__arg, label = "Enter value:", value = input[[write.fun__arg]]) 
    }
  })
  
  
  ### Data import:
  Dataset <- reactive({
    if (is.null(input$file)) {
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    args <- grep(paste0("^",input$readFunction,"__"), names(input), value = TRUE)
    
    argList <- list()
    for (i in seq_along(args))
    {
      argList[[i]] <- eval(parse(text=input[[args[i]]]))
    }
    names(argList) <- gsub(paste0("^",input$readFunction,"__"),"",args)
    
    argList <- argList[names(argList) %in% ArgNames()]
    
    Dataset <- as.data.frame(do.call(input$readFunction,c(list(input$file$datapath),argList)))
    return(Dataset)
  })
  
  # Select variables:
  output$varselect <- renderUI({
    
    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)
    
    # Variable selection:    
    selectInput("vars", "Variables to use:", names(Dataset()), multiple =TRUE)            
  })
  
  # Show table:
  output$table <- renderTable({
    
    if (is.null(input$vars) || length(input$vars)==0) return(NULL)
    
    Dataset2 = head(Dataset()[,input$vars,drop=FALSE], n = input$number.of.observations)
    
    return(Dataset2)
  })
  
  ### Download dump:
  
  output$downloadDump <- downloadHandler(
    filename = "Rdata.R",
    content = function(con) {
      
      assign(input$name, Dataset()[,input$vars,drop=FALSE])
      
      dump(input$name, con)
    }
  )
  
  ### Download save:
  
  output$downloadSave <- downloadHandler(
    filename = "Rdata.RData",
    content = function(con) {
      
      assign(input$name, Dataset()[,input$vars,drop=FALSE])
      
      save(list=input$name, file=con)
    }
  )
  
  
})

