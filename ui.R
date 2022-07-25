library(shiny)
library(shinyWidgets)

helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover({container: 'body'}); })")
      )
    ),
    tags$a(
      href = "#", class = "btn bt
      n-mini", `data-toggle` = "popover",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      HTML("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"16\" height=\"16\" fill=\"lightblue\" class=\"bi bi-question-circle-fill\" viewBox=\"0 0 16 16\">
   <path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.496 6.033h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286a.237.237 0 0 0 .241.247zm2.325 6.443c.61 0 1.029-.394 1.029-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94 0 .533.425.927 1.01.927z\"/>
 </svg>")
    )
  )
}

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
    # sliderInput(
    #   inputId = "K",
    #   label = "Disease prevalence:",
    #   min = 0.001,
    #   max = 0.2,
    #   step = 0.001,
    #   value = 0.01
    # )),
    sliderTextInput(
      inputId = "K",
      label = "Disease prevalence:",
      choices = unique(round(exp(seq(log(0.01), log(0.2), length = 500)), digits = 4)),
      grid = F, force_edges = T)),
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
        value = 0.3
      ))
    ),
  )),
  fluidRow(column(12, plotOutput(outputId = "distPlot", height = 600)))
)

slider_and_numeric <- function(id, label, min, max, step, value, helptext = "",
                               placement = "bottom") {
  if(length(step) > 1) {
    splitLayout(cellWidths = c("80%", "20%"),
      sliderTextInput(
        inputId = id,
        label = label,
        choices = step,
        grid = T, force_edges = T),
      helpPopup(NULL, helptext, placement = placement, c("hover")))
  }
  else {
    splitLayout(cellWidths = c("80%", "20%"),
                sliderInput(
                  inputId = id,
                  label = label,
                  min = min,
                  max = max,
                  step = step,
                  value = value
                ), 
                helpPopup(NULL, helptext, placement = placement, c("hover")))
  }
#            HTML(sprintf("<div title = \"%s\"> <svg xmlns=\"http://www.w3.org/2000/svg\" width=\"16\" height=\"16\" fill=\"currentColor\" class=\"bi bi-question-circle-fill\" viewBox=\"0 0 16 16\">
#   <path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.496 6.033h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286a.237.237 0 0 0 .241.247zm2.325 6.443c.61 0 1.029-.394 1.029-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94 0 .533.425.927 1.01.927z\"/>
# </svg></div>", helptext)))
           # ), HTML("<button type=\"button\" class=\"fa fa-question-circle\" data-toggle=\"popover\" title=\"Popover title\" data-content=\"And here's some amazing content. It's very engaging. Right?\">?</button>"))
}

calc_panel <- div(class = "well",
    #               tags$head(tags$script(HTML("$(function() {
    # setTimeout(function(){
    #   var vals = [0];
    #   var powStart = 0.01;
    #   var powStop = 0.2;
    #   while (powStart <= powStop) {
    #     vals.push(powStart);
    #     powStart *= 2;
    #   }
    #   $('#K2').data('ionRangeSlider').update({'values':vals})
    #   $('#K').data('ionRangeSlider').update({'values':vals})
    #   $('#K_1').data('ionRangeSlider').update({'values':vals})
    #   $('#K_2').data('ionRangeSlider').update({'values':vals})
    # }, 3)})"))),
    
    # tags$head(tags$script(HTML("$(function() {
    # setTimeout(function(){
    # $('#K2').data('ionRangeSlider').update({
    #        'prettify': function (num) { 
    #        return ((0.01*Math.pow(2, (num*7))).toFixed(2)); 
    #        }})
    # }, 5)})"))),
    tags$head(tags$script(("
    Number.prototype.mapLog = function (min, max) {
      const mapped = (this - min) * (Math.log(max) - Math.log(min)) / (max - min) + Math.log(min);
      return Math.exp(mapped);
    }
    $(function() {
    setTimeout(function(){
    
    # $('#K2').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #        return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
    #        }})

    $('#K_1').data('ionRangeSlider').update({
           'prettify': function (n) {
           return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
           }})
    $('#K_2').data('ionRangeSlider').update({
           'prettify': function (n) {
           return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
           }})
                           }, 2)})"))),
    
    # tags$head(tags$script(("
    # Number.prototype.mapLog = function (min, max) {
    #   const mapped = (this - min) * (Math.log(max) - Math.log(min)) / (max - min) + Math.log(min);
    #   return Math.exp(mapped);
    # }
    # (function() {
    # setTimeout(function(){
    # $('#K').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #         n.mapLog(this.min, this.max).toLocaleString('en-US');
    #        }})
    # $('#K2').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #         n.mapLog(this.min, this.max).toLocaleString('en-US');
    #        }})
    # $('#K_1').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #         n.mapLog(this.min, this.max).toLocaleString('en-US');
    #        }})
    # $('#K_2').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #         n.mapLog(this.min, this.max).toLocaleString('en-US');
    #        }})}, 2)})"))),
 
  fluidRow(column(4, slider_and_numeric("N2", "Number of embryos:", 2, 20, 1, 5, "The number of embryos available for selection."),
                  # slider_and_numeric("K2", "Disease prevalence:", 0.01, 1, 0.01, 0.5, NULL),
                  slider_and_numeric("K2", "Disease prevalence:", 0.01, 0.2, 
                                     sort(unique(c(seq(0.01, 0.2, 0.01), round(exp(seq(log(0.01), log(0.2), length = 500)), digits = 4)))), 0.01, "How prevalent is the disease in the population? 0.01 means 1% has the disease, while 0.2 means 20% of the population has the disease."),
                  slider_and_numeric("r2", "$$r^2:$$", 0, 1, NULL, 0.3, "The R-squared of the polygenic risk score model.")),
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
             slider_and_numeric("q2", "Quantile from which to exclude embryos:", 0, 1, 0.01, 0.3, paste("The quantile of the polygenic risk score from which we exclude embryos.")),
           ), fluidRow(column(8, offset = 2, htmlOutput("summary"), align = "center"))),
           column(4,
           conditionalPanel(
             condition = "input.type2 == 'Conditional'",
             slider_and_numeric("qf2", "Father's polygenic risk score quantile:", 0, 100, 1, 90, paste("Father's polygenic risk score quantile.")),
             slider_and_numeric("qm2", "Mother's polygenic risk score quantile:", 0, 100, 1, 90, paste("Mother's polygenic risk score quantile.")),
           ),
           conditionalPanel(
             condition = "input.type2 == 'Family History'",
             slider_and_numeric("h2", "$h^2:$", 0, 1, 0.01, 0.4, "The heritability of the disease."),
             checkboxInput("df2",
                           "Father has the disease"),
             checkboxInput("dm2",
                           "Mother has the disease"),
             slider_and_numeric("samples", "Number of monte carlo samples:", 5000, 300000, 10, 10000, "The number of simulations. Higher number will give a more accurate estimates, but might take longer to run.")))))

calc_two_traits <- div(class = "well", fluidRow(column(4,
                                                       slider_and_numeric("N_2", "Number of embryos:", 2, 20, 1, 5, "The number of embryos available for selection."),
                                                       slider_and_numeric("rho", '$\\rho$, the genetic correlation between the diseases:', -0.99, 0.99, 0.01, 0, "The genetic correlation between the two diseases."),
                                                       slider_and_numeric("samples_2", "Number of monte carlo samples:", 5000, 300000, 10, 10000, "The number of simulations. Higher number will give a more accurate estimates, but might take longer to run.")),
                                                column(4, 
                                                       slider_and_numeric("r2_1", "$$r^2 ~ \\text{disease 1:}$$", 0, 1, 0.01, 0.3, "The R-squared of the polygenic risk score model for the first disease."),
                                                       slider_and_numeric("r2_2", "$$r^2 ~ \\text{disease 2:}$$", 0, 1, 0.01, 0.3, "The R-squared of the polygenic risk score model for the second disease."),
                                                       fluidRow(column(8, offset = 2, htmlOutput("two_traits"), align = "center"))),
                                                column(4, 
                                                       slider_and_numeric("K_1", "Prevalence of disease 1:", 0.001, 0.2, unique(round(exp(seq(log(0.01), log(0.2), length = 500)), digits = 4)), 0.01, "How prevalent is the first disease in the population? 0.01 means 1% has the disease, while 0.2 means 20% of the population has the disease."),
                                                       slider_and_numeric("K_2", "Prevalence of disease 2:", 0.001, 0.2, unique(round(exp(seq(log(0.01), log(0.2), length = 500)), digits = 4)), 0.01, "How prevalent is the second disease in the population? 0.01 means 1% has the disease, while 0.2 means 20% of the population has the disease."))))



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


