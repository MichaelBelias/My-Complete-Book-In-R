
library(shiny)
ui <- fluidPage(
  sliderInput(inputId = "num",
              label = "Number of Random Normal Values",
              value = 550, min = 100, max = 1000),
  sliderInput(inputId = "breaks",
              label = "Number of Breaks in Histogram",
              value = 25, min = 10, max = 100),
  textInput(inputId = "title",
            label = "Write a title",
            value = "Histogram of Random Normal Values"),
  plotOutput(outputId = "hist"),
  verbatimTextOutput("stats")
)

server <- function(input, output) {
  data <- reactive({
    rnorm(input$num)
  })
  output$hist <- renderPlot({
    hist(data())
  })
  output$stats <- renderPrint({
    summary(data())
  })
}

shinyApp(server = server, ui = ui)