library(shiny)

about_panel <- verticalLayout(div(class = "well", h1("About", align = "center"),
                              HTML("<p>This application implements methods to calculate the expected risk reduction for a disease, when selecting embryos based on their polygenic risk scores(PRS). The calculations are based on the derivations in <a href = \"https://doi.org/10.7554/elife.64716\">Lencz etc.</a></p>"),
                              p("There are two selection strategies:"),
                              HTML("<ol>
                                   <li> Choosing the embryo with the lowest PRS out of all available embryos.</li>
                                   <li> Excluding embryos with PRS larger than some PRS quantile, and then picking a random embryo out of the remaining. In the case where there are no remaining embryos, we pick one at random.</li></ol>"),
                              h2("Plot"),
                              p("Plot allows you to generate graphs of the relative and absolute risk reduction as a function of a few variables, and depending on the embryo selection strategies."),
                              p("The variables are:"),
                              HTML("<ol>
                                   <li> $r^2$ - That is the R-squared of the PRS model.</li>
                                   <li> Disease prevalence - What is the prevalence (in percentage) of the disease in the population.</li>
                                   <li> Number of embryos - How many embryos are available to select from.</li>
                                   <li> Quantile from which to exclude - In the second strategy, we will select only embryos with PRS smaller than this specific quantile.</li>
                                   </ol>"),
                              h2("Calculator"),
                              p("Calculator is a more elaborate risk calculator based on the different variables, it also includes the ability to calculate the expected risk when conditioning on the parents' PRS or when conditioning on the parents' disease status. Note that in the latter case, it is based on simulation and might give a somewhat different answers between runs."),
                              p("When conditioning, we also have the following variables:"),
                              HTML("<ol>
                                   <li> Father's polygenic risk score quantile - pretty self explanatory, what is the father's quantile PRS score? </li>
                                   <li> Mother's polygenic risk score quantile - same as the previous one, but for the mother.</li>
                                   <li> $h^2$ - the heritability of the disease.</li>
                                   <li> Father has the disease - self explanatory.</li>
                                   <li> Mother has the disease - as the previous one.</li>
                                   <li> Number of monte carlo samples - the number of times the simulation is run. </li>
                                   </ol>"),
                              h2("Two diseases"),
                              p("Two diseases is a calculator of the risk of two diseases, given a certain correlation and when selecting the embryo with the lowest PRS for disease 1. It is based on a simulation as explained in the paper."),
                              p("The only new variable here is $\\rho$, the genetic correlation between the two diseases."),
                              h2("Reference"),
                              HTML("<ol>
                                   <li> <cite>Lencz, T., Backenroth, D., Granot-Hershkovitz, E., Green, A., Gettler, K., Cho, J. H., Weissbrod, O., Zuk, O., & Carmi, S. (2021). Utility of polygenic embryo screening for disease depends on the selection strategy. ELife, 10. <a href=\"https://doi.org/10.7554/elife.64716\">https://doi.org/10.7554/elife.64716</a> </cite> </li>
                                   </ol>")))

plot_panel <- div(class = "well",
  fluidRow(column(4,
    radioButtons("x_var", "Variable for x axis", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos"), 
                 choiceValues = c("r2", "Disease prevalence", "Number of embryos"),
                 selected = "r2",
                 inline = T),
    radioButtons(inputId = "lowestexclude",
      label = "Choose lowest risk embryo or exclude high risk embroys",
      choices = c("Lowest", "Exclude"), inline = T
    )),
    column(4, tags$div(title = "The R-squared of the polgenic risk score model.", sliderInput(
      inputId = "r",
      label = "$$r^2:$$",
      min = 0,
      max = 1,
      step = 0.01,
      value = 0.3
    )),
    sliderInput(
      inputId = "K",
      label = "Disease prevalence:",
      min = 0.01,
      max = 1,
      step = 0.01,
      value = 0.05
    )),
    column(4, sliderInput(
      inputId = "N",
      label = "Number of embryos:",
      min = 2,
      max = 20,
      step = 1,
      value = 5
    ),
    conditionalPanel(
      condition = "input.lowestexclude == 'Exclude'",
      tags$div(title = "The quantile of the polygenic risk score from which we exclude embryos.", sliderInput(
        inputId = "q",
        label = "Quantile from which to exclude:",
        min = 0,
        max = 1,
        step = 0.01,
        value = 0.5
      ))
    ),
  )),
  fluidRow(column(12, plotOutput(outputId = "distPlot", height = 600)))
)

slider_and_numeric <- function(id, label, min, max, step, value, helptext = "") {
  column(width = 12,
         tags$div(title = helptext,
                  sliderInput(
                    inputId = id,
                    label = label,
                    min = min,
                    max = max,
                    step = step,
                    value = value
                  )),
         numericInput(paste0("input_", id), label = "", value = value, min = min, max = max, step = step), offset = 0)
}

calc_panel <- div(class = "well",
  fluidRow(column(4, slider_and_numeric("N2", "Number of embryos:", 2, 1000, 1, 5, NULL),
                  slider_and_numeric("K2", "Disease prevalence:", 0.01, 1, 0.01, 0.5, NULL),
                  slider_and_numeric("r2", "$$r^2:$$", 0, 1, 0.01, 0.5, "The R-squared of the polygenic risk score model.")),
           column(4, radioButtons(
             inputId = "lowestexclude2",
             label = "Choose lowest risk embryo or exclude high risk embroys",
             choices = c("Lowest", "Exclude"), inline = T),
           radioButtons(
             inputId = "type2",
             label = "Should we condition on the parents information?",
             # choices = c("Risk reduction", "Conditional", "Family History"), inline = T
             choiceValues = c("Risk reduction", "Conditional", "Family History"),
             choiceNames = c("No conditioning", "Conditional on parents polygenic risk score", "Conditional on parents disease status"), 
           ),
           conditionalPanel(
             condition = "input.lowestexclude2 == 'Exclude'",
             slider_and_numeric("q2", "Quantile from which to exclude embryos:", 0, 1, 0.01, 0.5, paste("The quantile of the polygenic risk score from which we exclude embryos.")),
           ), fluidRow(column(8, offset = 2, htmlOutput("summary"), align = "center"))),
           column(4,
           conditionalPanel(
             condition = "input.type2 == 'Conditional'",
             slider_and_numeric("qf2", "Father's polygenic risk score quantile:", 0, 100, 1, 5, paste("Father's PRS quantile.")),
             slider_and_numeric("qm2", "Mother's polygenic risk score quantile:", 0, 100, 1, 5, paste("Mother's PRS quantile.")),
           ),
           conditionalPanel(
             condition = "input.type2 == 'Family History'",
             slider_and_numeric("h2", "$h^2:$", 0, 1, 0.01, 0.5, "The heritability of the disease."),
             checkboxInput("df2",
                           "Father has the disease"),
             checkboxInput("dm2",
                           "Mother has the disease"),
             slider_and_numeric("samples", "Number of monte carlo samples:", 5000, 300000, 10, 10000, paste("The number of simulations"))))))

calc_two_traits <- div(class = "well", fluidRow(column(4,
                                                       slider_and_numeric("N_2", "Number of embryos:", 2, 100, 1, 5, NULL),
                                                       slider_and_numeric("rho", '$\\rho$, the genetic correlation between the diseases:', -0.99, 0.99, 0.01, 0, "The correlation between the two diseases."),
                                                       slider_and_numeric("samples_2", "Number of monte carlo samples:", 5000, 300000, 10, 10000, "The number of simulations.")),
                                                column(4, 
                                                       slider_and_numeric("r2_1", "$$r^2 ~ \\text{disease 1:}$$", 0, 1, 0.01, 0.5, "The R-squared of the polygenic risk score model for the first disease"),
                                                       slider_and_numeric("r2_2", "$$r^2 ~ \\text{disease 2:}$$", 0, 1, 0.01, 0.5, "The R-squared of the polygenic risk score model for the second disease"),
                                                       fluidRow(column(8, offset = 2, htmlOutput("two_traits"), align = "center"))),
                                                column(4, 
                                                       slider_and_numeric("K_1", "Prevalence of disease 1:", 0, 1, 0.01, 0.5, NULL),
                                                       slider_and_numeric("K_2", "Prevalence of disease 2:", 0, 1, 0.01, 0.5, NULL))))



ui <- fluidPage(
  tags$head(HTML("<script type=\"text/x-mathjax-config\">
  MathJax.Hub.Config({
    tex2jax: {
      inlineMath: [ ['$','$'] ],
      processEscapes: true
    }
  });
</script>")),
  shinyjs::useShinyjs(),
  withMathJax(),
  tabsetPanel(
    type = "pills",
    tabPanel("About", about_panel),
    tabPanel("Plot", plot_panel),
    tabPanel("Calculator", calc_panel),
    tabPanel("Two diseases", calc_two_traits)
  ))


