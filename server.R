library(shiny)
library(ggplot2)

source("EmbryoSelection.R")
server <- function(input, output, session) {
  accepeted <- F
  modal <- modalDialog(h1("Important notes:"), 
                       HTML("<ol>
       <li>Our risk estimates are based on the liability threshold model from statistical genetics. They are not directly based on real epidemiological data.</li>
       <li>Under the model, a disease is assumed to have an underlying, continuous liability. The liability is assumed to be normally distributed and to represent the sum of additive genetic and environmental risk factors. Under the model, an individual is affected if his or her liability exceeds a threshold.</li>
       <li>The model assumes that each adult individual is either affected or unaffected. Both the prevalence parameter and the risk estimated by the application represent the proportion of affected adults.</li>
       <li>The model assumes that a PRS explains a proportion R^2 of the variance in the liability in the population. This parameter, which quantifies the accuracy of the score, should already incorporate any reduction in the accuracy of the score due to genotyping errors, an application of the score in populations of non-European descent, elimination of environmental correlates of the score in the setting of comparing sibling embryos, and other factors.</li>
       <li>All embryos in the model are assumed to have the potential to lead to live birth when transferred. If you instead have, for example, the number of blastocysts, please adjust it accordingly based on specifics of your given case. We assume a single embryo transfer.</li>
       <li>We assume that the parents screen the embryos for a <b>single disease</b>. We assume random mating (i.e., no assortative mating related to the disease being screened). The calculator does not provide information on the likelihood of increasing the risk of other diseases due to selection. See our paper for more details.</li>
       <li>Finally, we note that screening IVF embryos with polygenic risk scores may be associated with ethical, social, and legal problems. Please see our paper for more details.</li>
       </ol>"),
                       p("Please confirm that:"),
                       HTML("<ol>
       <li>You read and understood the above notes.</li>
       <li>You understand that the application is intended for research purposes only and is not intended to guide clinical decision making.</li>
       </ol>"), checkboxInput("accept", "Accept"), easyClose = F, footer = NULL, size = "l")

  observeEvent(input$accept, {
    req(input$accept)
    removeModal()
    accepeted <<- T
  })

  observeEvent(input$lowestexclude, {
    if (input$lowestexclude == "Lowest") {
      # withMathJax(updateSelectInput(session, "x_var", choices = c("r2", "Disease prevalence", "Number of embryos")))
      withMathJax(updateRadioButtons(session, "x_var", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos"), 
                                     choiceValues = c("r2", "Disease prevalence", "Number of embryos"), inline = T))
    }
    else {
      # withMathJax(updateSelectInput(session, "x_var", choices = c("r2", "Disease prevalence", "Number of embryos", "Quantile")))
      withMathJax(updateRadioButtons(session, "x_var", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos", "Percentile"),
                                     choiceValues = c("r2", "Disease prevalence", "Number of embryos", "Percentile"), inline = T))
    }
  })
  
  observeEvent(input$x_var, {
    if (input$x_var == "Number of embryos") {
      shinyjs::hide("N")
      shinyjs::show("r")
      shinyjs::show("K")
      shinyjs::show("q")
    }
    else if (input$x_var == "Disease prevalence") {
      shinyjs::hide("K")
      shinyjs::show("r")
      shinyjs::show("N")
      shinyjs::show("q")
    }
    else if (input$x_var == "r2") {
      shinyjs::hide("r")
      shinyjs::show("N")
      shinyjs::show("K")
      shinyjs::show("q")
    }
    else if (input$x_var == "Percentile") {
      shinyjs::hide("q")
      shinyjs::show("N")
      shinyjs::show("K")
      shinyjs::show("r")
    }
  })
  
  filter_input <- debounce(reactive({
    list(r2 = input$r,
         K2 = input$K,
         N = input$N,
         q = input$q)
  }), 500)

  output$distPlot <- renderPlot({
    if (!accepeted) showModal(modal)
    
    filts <- filter_input()
    # output$distPlot <- renderPlotly({
    # selectInput("x_var", "Variable for x axis", choices = c("r2", "K", "N")),
    if (input$lowestexclude == "Lowest") {
      subtitle <- "Lowest PRS strategy"
    }
    else {
      subtitle <- "Exclude high PRS strategy"
    }
    selected_x <- input$x_var
    if(selected_x == "Number of embryos") {
      x <- 2:20
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(filts$r, filts$K, n = x))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(r2 = filts$r, K = filts$K, q = filts$q, n = x))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      
      # subtitle <- "Lowest strategy"
      x_lab <- "Number of embryos"
    }
    else if(selected_x == "Disease prevalence") {
      x <- exp(seq(log(0.001), log(0.3), length = 50))
      # y <- sapply(x, function(x) risk_reduction_lowest(input$r, x, n = input$N))
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(filts$r, x, n = filts$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * x
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(r2 = filts$r, K = x, q = filts$q, n = filts$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * x
      }
      # subtitle <- "Lowest strategy"
      x_lab <- "Disease prevalence"
    }
    else if(selected_x == "r2" || selected_x == "$$R^2$$") {
      x <- seq(0.01, 1, length = 50)
      # y <- sapply(x, function(x) risk_reduction_lowest(x, input$K, n = input$N))
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(x, filts$K, n = input$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(x, filts$K, filts$q, n = filts$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      # subtitle <- "Lowest strategy"
      # x_lab <- "PRS r^2"
      x_lab <- expression(R^2)
    }
    else if (selected_x == "Percentile") {
      # q
      x <- seq(0.01, 1, length = 50)
      y <- sapply(x, function(x) risk_reduction_exclude(filts$r, filts$K, x, n = filts$N))
      # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      # subtitle <- "Exclude strategy"
      x_lab <- "Percentile to exclude"
    }
    else {
      print("Error!")
      print(selected_x)
    }
    
    if(selected_x == "Disease prevalence") {
      ggplot() +
        geom_point(mapping = aes(x, y, color = "Relative risk"), size = 3) +
        geom_line(mapping = aes(x, y, color = "Relative risk")) +
        geom_point(mapping = aes(x, y*x, color = "Absolute risk"), size = 3) +
        geom_line(mapping = aes(x, y*x, color = "Absolute risk")) +
        theme_minimal() +
        labs(
          title = "Risk reduction",
          subtitle = subtitle,
          x = x_lab,
          y = "Risk reduction"
        ) +
        theme(
          plot.subtitle = element_text(size = 22, hjust = 0.5),
          plot.title = element_text(size = 25, hjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          panel.background = element_rect(fill = "#f5f5f5", color = NA),
          plot.background = element_rect(fill = "#f5f5f5", color = NA),
          panel.grid.minor = element_line(colour = "darkgray"),
          panel.grid.major = element_line(colour = "darkgray"),
          legend.text = element_text(size = 16),
          legend.title = element_blank()
        ) +
        scale_y_continuous(name = "Risk reduction")
    }
    else {
      ggplot(mapping = aes(x, y)) +
        geom_point(size = 3) +
        geom_line() +
        theme_minimal() +
        labs(
          title = "Risk reduction",
          subtitle = subtitle,
          x = x_lab,
          y = "Risk reduction"
        ) +
        theme(
          plot.subtitle = element_text(size = 22, hjust = 0.5),
          plot.title = element_text(size = 25, hjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          panel.background = element_rect(fill = "#f5f5f5", color = NA),
          plot.background = element_rect(fill = "#f5f5f5", color = NA),
          panel.grid.minor = element_line(colour = "darkgray"),
          panel.grid.major = element_line(colour = "darkgray")
        ) +
        scale_y_continuous(name = "Relative risk reduction",
                           sec.axis = sec_axis(~.*ifelse(selected_x == "K", x, input$K), name = "Absolute risk reduction")) 
    }
  })
  
  output$summary <- renderPrint({
    if (!accepeted) showModal(modal)
    # print(input$K2)
    cat("<div class = \"alert alert-info\">")
    if (input$type2 == "Risk reduction") {
      # temp <- risk_reduction_lowest(input$r2, K = input$K2, n = input$N2)
      temp <- risk_reduction_lowest(input$r2, K = input$K2, 
                                    n = input$N2)
      if (input$lowestexclude2 != "Lowest") {
        # temp <- risk_reduction_exclude(input$r2,
        #                          K = input$K2,
        #                          q = input$q2,
        #                          n = input$N2)
        temp <- risk_reduction_exclude(input$r2,
                                       K = input$K2,
                                       q = input$q2,
                                       n = input$N2)
      }
      cat(
        sprintf(
          "<p><u>Baseline risk</u>: <strong>%.4f</strong>\n</p>",
          input$K2
        )
      )
      cat(
        sprintf(
          "<p><u>Risk for specific strategy</u>: <strong>%.4f</strong>\n</p>",
           input$K2 - temp * input$K2
        )
      )
      cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.4f</strong>\n</p>", temp))
      cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.4f</strong></p>", temp * input$K2))
      
      cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f</strong></p>", ceiling(1/(temp * input$K2))))
    }
    else if (input$type2 == "Conditional") {
      temp <-
        risk_reduction_lowest_conditional(
          input$r2,
          K = input$K2,
          n = input$N2,
          qf = 1-input$qf2,
          qm = 1-input$qm2,
          relative = T,
          parental_avg_given = F
        )
      if (input$lowestexclude2 != "Lowest") {
        temp <-
          risk_reduction_exclude_conditional(
            input$r2,
            K = input$K2,
            q = input$q2,
            n = input$N2,
            qf = 1-input$qf2,
            qm = 1-input$qm2,
            relative = T
          )
      }
      
      cat(
        sprintf(
          "<p><u>Baseline risk</u>: <strong>%.4f</strong>\n</p>",
          temp$baseline
        )
      )
      cat(
        sprintf(
          "<p><u>Risk for specific strategy</u>: <strong>%.4f</strong>\n</p>",
          temp$risk
        )
      )
      cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.4f</strong>\n</p>", temp$rr))
      cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.4f</strong></p>", temp$rr * temp$baseline))
      cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f</strong></p>", ceiling(1/(temp$rr * temp$baseline))))
    }
    else {
      if (input$r2 >= input$h2) {
        updateNumericInput(session, "input_r2", value = input$h2)
        updateSliderInput(session, "r2", value = input$h2)
      }
      
      
      if (input$r2 <= input$h2) {
        if (input$lowestexclude2 != "Lowest") {
          temp <-
            risk_reduction_exclude_family_history(
              input$r2,
              h2 = input$h2,
              K = input$K2,
              q = input$q2,
              n = input$N2,
              input$df2,
              input$dm2,
              n_samples = input$samples
            )
        }
        else {
          temp <-
            risk_reduction_lowest_family_history(
              input$r2,
              h2 = input$h2,
              K = input$K2,
              n = input$N2,
              input$df2,
              input$dm2,
              n_samples = input$samples
            )
        }
        cat(
          sprintf(
            "<p><u>Baseline risk</u>: <strong>%.4f</strong>\n</p>",
            temp[1]
          )
        )
        cat(
          sprintf(
            "<p><u>Risk for specific strategy</u>: <strong>%.4f (%.4f)</strong>\n</p>",
            temp[2],
            temp[5]
          )
        )
        
        cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.4f (%.4f) </strong>\n</p>", temp[3], temp[5] / temp[1]))
        cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.4f (%.4f)</strong>\n</p>", temp[4], temp[5]))
        
        cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f (%.4f)</strong>\n</p>", ceiling(1/temp[4]), temp[5] / temp[4]^2))
        
        cat(sprintf("<p style = \"color:red\">Based on %d draws. Estimated standard deviation in parentheses.</p>", input$samples))
        if(temp[1] < temp[2]) {
          cat(sprintf("<b><p style = \"color:red\">Baseline risk is smaller than the strategy risk. Either the sample size is too small, or the baseline and strategy risk are almost identical.</p></b>"))
        }
        # cat("Note that the standard deviation of the relative risk reduction can be much higher, and a larger sample size should be used to estimate it well.")
      }
      else {
        updateNumericInput(session, "input_r2", value = input$h2)
        updateSliderInput(session, "r2", value = input$h2)
        # withMathJax(cat("$$r^2\text{ can't be higher than }h^2$$"))
      }
    }
    cat("</div>")
  })
  
  output$two_traits <- renderPrint({
    cat("<div class = \"alert alert-info\">")
    temp <-
      simulate_lowest_risk_two_traits(input$r2_1, input$r2_2, input$rho, input$K_1, input$K_2, input$N_2, input$samples_2)
    cat(sprintf("<p><u>Relative risk reduction for disease 1</u>: <strong>%.4f</strong>\n</p>", temp[1]))
    cat(sprintf("<p><u>Absolute risk reduction for disease 1</u>: <strong>%.4f</strong>\n</p>", temp[2]))
    cat(sprintf("<p><u>Couples needed to screen for disease 1</u>: <strong>%.0f</strong>\n</p>", ceiling(1/temp[2])))
    
    cat(sprintf("<p><u>Relative risk reduction for disease 2</u>: <strong>%.4f</strong>\n</p>", ifelse(input$rho == 0, 0, temp[3])))
    cat(sprintf("<p><u>Absolute risk reduction for disease 2</u>: <strong>%.4f</strong>\n</p>", ifelse(input$rho == 0, 0, temp[4])))
    cat(sprintf("<p><u>Couples needed to screen for disease 2</u>: <strong>%.0f</strong>\n</p>", ifelse(input$rho == 0, 0, ceiling(1/temp[4]))))
    
    cat(sprintf("<p style = \"color:red\">Based on %d simulations.</p>", input$samples_2))
    cat("</div>")
  })
}

# shinyApp(ui = ui, server = server)
# runApp(launch.browser = F)