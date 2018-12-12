
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Simulations"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId="beta.0",label = "Beta 0",value = -2,min = -Inf,max = Inf),
      numericInput(inputId="beta.1",label = "Beta 1",value = 5,min = -Inf,max = Inf),
      numericInput(inputId="beta.2",label = "Beta 2",value = 20,min = -Inf,max = Inf),
      numericInput(inputId="n",label = "Number of observations",value = 100,min = 1,max = Inf,step = 1),
      actionButton("go","Go")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
