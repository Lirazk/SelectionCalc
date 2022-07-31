# Done
# Change r^2 to R^2 to be consistent (in the plot axis)
# Number of embryos - Change max to 10
# Remove two diseases
# Number needed to screen (change to number of couples) - use ceiling
# Monte carlo question mark - a more accurate estimate (without s)
# Plot - stronger grid lines
# Plot - larger points, quantiles from which to exclude should be > 0
# Plot - the quantile from which to exclude, add example, "PRS above", add question mark. "if all embryos are high risk, a random embryo is chosen"
# Prevalence - start from 1/1000, no ticks, to 0.3
# Also change box site so every number would fit
# Change quantile to percentile
# Add space before parentheses
# Plot - Prevalence by RR, does it make sense? shouldn't the abs risk and RR be different?
# Add question marks when plot
# Disease prevalence question mark - means that 1% of the population have the disease
# The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the square. Typically equal 0.05-0.1.
# h^2 - add "only relevant when conditioning on the parents' diseases status"
# parents' polygenic risk score 
# "Based on x simulations."
# Explain father/mother PRS, and change it so 1% is the top 1% risk

# TODO
# Call it: PGT-P outcome calculator (single disease).
# Add tests.

library(shiny)
library(shinyWidgets)

slider_and_numeric <- function(id, label, min, max, step, value, helptext = "",
                               placement = "bottom") {
  if(length(step) > 1) {
    div(id = id, 
        splitLayout(cellWidths = c("80%", "20%"),
                    sliderTextInput(
                      inputId = id,
                      label = label,
                      choices = step,
                      # grid = T, force_edges = T),
                      grid = F, force_edges = T),
                    helpPopup(NULL, helptext, placement = placement, c("hover"))))
    # splitLayout(cellWidths = c("80%", "20%"),
    #             sliderTextInput(
    #               inputId = id,
    #               label = label,
    #               choices = step,
    #               # grid = T, force_edges = T),
    #               grid = F, force_edges = T),
    #             helpPopup(NULL, helptext, id = id, placement = placement, c("hover")))
  }
  else {
    div(id = id, 
    splitLayout(cellWidths = c("80%", "20%"),
                sliderInput(
                  inputId = id,
                  label = label,
                  min = min,
                  max = max,
                  step = step,
                  value = value
                ), 
                helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
  #            HTML(sprintf("<div title = \"%s\"> <svg xmlns=\"http://www.w3.org/2000/svg\" width=\"16\" height=\"16\" fill=\"currentColor\" class=\"bi bi-question-circle-fill\" viewBox=\"0 0 16 16\">
  #   <path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.496 6.033h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286a.237.237 0 0 0 .241.247zm2.325 6.443c.61 0 1.029-.394 1.029-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94 0 .533.425.927 1.01.927z\"/>
  # </svg></div>", helptext)))
  # ), HTML("<button type=\"button\" class=\"fa fa-question-circle\" data-toggle=\"popover\" title=\"Popover title\" data-content=\"And here's some amazing content. It's very engaging. Right?\">?</button>"))
}


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
                              HTML("<p>This application allows users to predict the expected risk reduction when selecting an IVF embryo for transfer based on polygenic risk scores (PRS) for a single disease. The model is described in <a href = \"https://doi.org/10.7554/elife.64716\" target=\"_blank\" rel=\"noopener noreferrer\">Lencz et al.</a></p>"),
                              p("We provide estimates under two selection strategies."),
                              HTML("<ol>
                                   <li>Lowest risk prioritization: Selecting the embryo with the lowest PRS among all available embryos.</li>
                                   <li>High-risk exclusion: Excluding embryos with a PRS above a “high-risk” cutoff, and then selecting an embryo at random among the remaining embryos. In the case all embryos are high-risk, a random embryo is selected.</li></ol>"),
                              h2("Plot mode"),
                              p("Use the Plot mode to generate graphs of the relative and absolute risk reductions vs the parameters of the problem for each of the two selection strategies."),
                              p("The parameters are:"),
                              HTML("<ol>
                                   <li>$R^2$: The proportion of variance in liability to the disease explained by the PRS, which is a measure of the accuracy of the score. $R^2$ is currently between 5-10% for most common polygenic diseases.</li>
                                   <li>The disease prevalence: The proportion of individuals in the (adult) population affected by the disease.</li>
                                   <li>The number of embryos: The number of viable embryos from which a single embryo is selected for transfer. Important note: the model assumes that each embryo transferred will be born. This parameter should therefore correspond to the number of live births expected from the given cycle.</li>
                                   <li>Quantile from which to exclude: In the “high-risk exclusion” strategy, this is the cutoff that defines embryos as high-risk. Embryos with PRS above that quantile are excluded. For example, if the parameter equals 10%, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded.</li>
                                   </ol>"),
                              h2("Calculator mode"),
                              p("Use the calculator mode to present the baseline risk and the risk reduction given a set of parameters. The calculator outputs the risk with and without embryo selection, the relative and the absolute risk reduction, and the number of couples needed to screen their embryos to prevent a single future case. The calculator also allows users to compute risk estimates when conditioning on the parental PRS percentiles or disease status."),
                              p("When conditioning on the parental PRS percentiles, the following parameters are required:"),
                              HTML("<ol>
                                   <li>Father's PRS percentile: For example, if this parameter equal 0.05, the PRS of the father is at the top 5% of the distribution of PRS in the population.</li>
                                   <li>Mother's PRS percentile: An analogous parameter for the mother.</li>
                                   </ol>"),
                              p("When conditioning on the parental disease status, the following parameters are required:"),
                              HTML("<ol>
                                   <li>Father has the disease: yes/no</li>
                                   <li>Mother has the disease: yes/no</li>
                                   <li>$h^2$: The heritability of the disease, i.e., the proportion of variance in the liability of the disease explained by additive genetic factors.</li>
                                   <li>Number of Monte-Carlo draws: when conditioning on the parental disease status, the risk is estimated using a Monte-Carlo method. With more draws, the risk estimate becomes more accurate but is slower to compute.</li>
                                   </ol>"),
                              # h2("Two diseases"),
                              # p("Two diseases is a calculator of the risk of two diseases, given a certain correlation and when selecting the embryo with the lowest PRS for disease 1. It is based on a simulation as explained in the paper."),
                              # p("The only new variable here is $\\rho$, the genetic correlation between the two diseases."),
                              h2("Reference"),
                              HTML("<ol>
                                   <li> <cite>Lencz, T., Backenroth, D., Granot-Hershkovitz, E., Green, A., Gettler, K., Cho, J. H., Weissbrod, O., Zuk, O., & Carmi, S. (2021). Utility of polygenic embryo screening for disease depends on the selection strategy. eLife, 10. <a href=\"https://doi.org/10.7554/elife.64716\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.7554/elife.64716</a> </cite> </li>
                                   </ol>"),
                              h2("Contact"),
                              p("Please contact us if you find an error or have any suggestion."),
                              HTML("<p>Shai Carmi, <a href=\"mailto: shai.carmi@huji.ac.il\">shai.carmi@huji.ac.il</a></p>"),
                              HTML("<p>Liraz Klausner, <a href=\"mailto: liraz.klausner@mail.huji.ac.il\">liraz.klausner@mail.huji.ac.il</a></p>"),
                              p("Braun School of Public Health, The Hebrew University of Jerusalem")),
                              p(paste("Last update", max(format(max(file.info("ui.R")$mtime,
                                                                    file.info("server.R")$mtime,
                                                                    file.info("EmbryoSelection.R")$mtime), "%d-%m-%Y")))))

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
    # column(4, tags$div(title = "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1.", sliderInput(
    #   inputId = "r",
    #   label = "$$R^2:$$",
    #   min = 0,
    #   max = 1,
    #   step = 0.01,
    #   value = 0.3
    # )),
    column(4, slider_and_numeric("r", "$$R^2:$$", 0.01, 1, NULL, 0.05, 
                                 "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
    # sliderInput(
    #   inputId = "K",
    #   label = "Disease prevalence:",
    #   min = 0.001,
    #   max = 0.2,
    #   step = 0.001,
    #   value = 0.01
    # )),
    
    # sliderTextInput(
    #   inputId = "K",
    #   label = "Disease prevalence:",
    #   choices = unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)),
    #   grid = F, force_edges = T)),
  
    slider_and_numeric("K", "Disease prevalence:", 0.001, 0.3, 
                       sort(unique(c(seq(0.001, 0.3, 0.001), round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 0.001, 
                       "How prevalent is the disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."),
    # column(4, sliderInput(
    #   inputId = "N",
    #   label = "Number of embryos:",
    #   min = 2,
    #   max = 10,
    #   step = 1,
    #   value = 5
    # ),
    slider_and_numeric("N", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection."),
    conditionalPanel(
      condition = "input.lowestexclude == 'Exclude'",
      # tags$div(title = "The percentile of the polygenic risk score from which we exclude embryos.", sliderInput(
      #   inputId = "q",
      #   label = "Percentile from which to exclude:",
      #   min = 0.01,
      #   max = 1,
      #   step = 0.01,
      #   value = 0.3
      # ))
      slider_and_numeric("q", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded. If no embryo is below the threshold we select one randomly."))
    )),
  ),
  fluidRow(column(12, plotOutput(outputId = "distPlot", height = 600)))
)

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
 
  fluidRow(column(4, slider_and_numeric("N2", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection."),
                  # slider_and_numeric("K2", "Disease prevalence:", 0.01, 1, 0.01, 0.5, NULL),
                  slider_and_numeric("K2", "Disease prevalence:", 0.001, 0.3, 
                                     sort(unique(c(seq(0.001, 0.3, 0.001), round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 0.001, "How prevalent is the disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."),
                  slider_and_numeric("r2", "$$R^2:$$", 0.01, 1, NULL, 0.05, "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1.")),
           column(4, radioButtons(
             inputId = "lowestexclude2",
             label = "Choose lowest risk embryo or exclude high risk embroys",
             choices = c("Lowest", "Exclude"), inline = T),
           radioButtons(
             inputId = "type2",
             label = "Should we condition on the parents' information?",
             # choices = c("Risk reduction", "Conditional", "Family History"), inline = T
             choiceValues = c("Risk reduction", "Conditional", "Family History"),
             choiceNames = c("No conditioning", "Conditional on the parents' polygenic risk score", "Conditional on the parents' disease status"), 
           ),
           conditionalPanel(
             condition = "input.lowestexclude2 == 'Exclude'",
             slider_and_numeric("q2", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of the PRS in the population will be excluded. If no embryo is below the threshold we select one randomly.")),
           ), fluidRow(column(10, offset = 2, htmlOutput("summary"), align = "center"))),
           column(4,
           conditionalPanel(
             condition = "input.type2 == 'Conditional'",
             slider_and_numeric("qf2", "Father's polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.05, the PRS of the father is at the top 5% of the distribution of the PRS in the population.")),
             slider_and_numeric("qm2", "Mother's polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.05, the PRS of the mother is at the top 5% of the distribution of the PRS in the population.")),
           ),
           conditionalPanel(
             condition = "input.type2 == 'Family History'",
             slider_and_numeric("h2", "$h^2:$", 0.01, 1, 0.01, 0.4, "The heritability of the disease. Only relevant when conditioning on the parents' disease status"),
             checkboxInput("df2",
                           "Father has the disease"),
             checkboxInput("dm2",
                           "Mother has the disease"),
             slider_and_numeric("samples", "Number of monte carlo draws:", 5000, 300000, 10, 10000, "The number of simulations. Higher number will give a more accurate estimate, but might take longer to run.")))))

calc_two_traits <- div(class = "well", fluidRow(column(4,
                                                       slider_and_numeric("N_2", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection."),
                                                       slider_and_numeric("rho", '$\\rho$, the genetic correlation between the diseases:', -0.99, 0.99, 0.01, 0, "The genetic correlation between the two diseases."),
                                                       slider_and_numeric("samples_2", "Number of monte carlo draws:", 5000, 300000, 10, 10000, "The number of simulations. Higher number will give a more accurate estimate, but might take longer to run.")),
                                                column(4, 
                                                       slider_and_numeric("r2_1", "$$R^2 ~ \\text{disease 1:}$$", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the first disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       slider_and_numeric("r2_2", "$$R^2 ~ \\text{disease 2:}$$", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the second disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       fluidRow(column(8, offset = 2, htmlOutput("two_traits"), align = "center"))),
                                                column(4, 
                                                       slider_and_numeric("K_1", "Prevalence of disease 1:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the first disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."),
                                                       slider_and_numeric("K_2", "Prevalence of disease 2:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the second disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."))))



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
    # tabPanel("Two diseases", calc_two_traits)
  ))


