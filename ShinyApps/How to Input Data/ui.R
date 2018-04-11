
shinyUI(pageWithSidebar(
  
  # Header:
  headerPanel("R data reader"),
  
  # Input in sidepanel:
  sidebarPanel(
    tags$style(type='text/css', ".well { max-width: 20em; }"),
    # Tags:
    tags$head(
      tags$style(type="text/css", "select[multiple] { width: 100%; height:10em}"),
      tags$style(type="text/css", "select { width: 100%}"),
      tags$style(type="text/css", "input { width: 19em; max-width:100%}")
    ),
    
    # Select filetype:
    selectInput("readFunction", "Function to read data:", c(
      # Base R:
      "read_sav",
      "read.table",
      "read.csv",
      "read.csv2",
      "read.delim",
      "read.delim2",
      
      # foreign functions:
      "read.arff",
      "read.dbf",
      "read.dta",
      "read.epiiinfo",
      "read.mtp",
      "read.octave",
      "read.ssd",
      "read.systat",
      "read.xport",
      
      # Advanced functions:
      "read.xlsx",
      "read.xlsx2",
      "scan",
      "readLines"
    )),
    
    # Argument selecter:
    htmlOutput("ArgSelect"),
    # Argument field:
    htmlOutput("ArgText"),
    
    
    # Upload data:
    fileInput("file", "Upload data-file:"),
    
    # Variable selection:
    htmlOutput("varselect"),
    
    textInput("name","Dataset name:","Data"),
    
    br(),
    
    # Select filetype:
    selectInput("write.Function", "Function to read data:", c(
      # Base R:
      "write_sav",
      "write.table",
      "write.csv",
      "write.csv2",
      "write.delim",
      "write.delim2",
      
      # foreign functions:
      "write.arff",
      "write.dbf",
      "write.dta",
      "write.epiiinfo",
      "write.mtp",
      "write.octave",
      "write.ssd",
      "write.systat",
      "write.xport",
      
      # Advanced functions:
      "write.xlsx",
      "write.xlsx2",
      "scan",
      "readLines"
    )),
    # write.Argument selecter:
    htmlOutput("write.ArgSelect"),
    # write.Argument field:
    htmlOutput("write.ArgText"),
    
    downloadLink('downloadDump', 'Download source'),
    downloadLink('downloadSave', 'Download binary'),
    sliderInput("number.of.observations", "observations to show" ,min = 1 , max=100, value = 10)
    
  ),
  
  # Main:
  mainPanel(
    tableOutput("table")
  )
))
