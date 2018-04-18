fluidPage(
  titlePanel("Uploading Files"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose file to upload',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv',
                  '.xlsx',
                  '.xlsm'
                )
      ),
      tags$hr(),
      p()
    ),
    mainPanel(
      
    )
  )
,
  # App title ----
  titlePanel("Choose your combination"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      
      # Input: Sensitivity of PLND with step value  0.05----
      sliderInput("PLND_sens", "Sensitivity of PLND:",
                  min = 0, max = 1,
                  value = 0.5, step = 0.05),
      
      # Input: Sensitivity of the new test with step value  0.05----
      sliderInput("New_Sens", "Sensitivity of the new test:",
                  min = 0, max = 1,
                  value = 0.5, step = 0.05),
      
      # Input: Sensitivity of the new test with step value  0.05----
      sliderInput("New_Spes", "Specificity of the new test:",
                  min = 0, max = 1,
                  value = 0.5, step = 0.05)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Table summarizing the values entered ----
      tableOutput("values"),
      tableOutput('Dataset'),
      tableOutput('table'),
      plotOutput("plot1"),
      plotOutput("plot2")
      
      
    )
  )
)

