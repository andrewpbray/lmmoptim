library(shiny)

fluidPage(fluidRow(plotOutput('bbplot', height = "400px")),
          fluidRow(style = "padding-bottom: 20px;",
                   column(2, numericInput("eps", label = "epsilon", value = .1, min = 0, max = .5, step = .1)),
                   column(2, numericInput("maxiter", label = "iterations", min = 1, max = 10, value = 8,
                                         step = 1)),
                   column(4, sliderInput("it", label = "stage", min = 1, max = 24, value = 1,
                                         step = 1, animate = animationOptions(loop = F, interval=2300))),
                   column(4, radioButtons("radio", label = "",
                                          choices = list("branch" = "option1", "bound" = "option2", "kill" = "option3"),
                                          selected = "option2"))
          )
)
