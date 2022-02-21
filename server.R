library(shiny)
library(ggplot2)

source("EmbryoSelection.R")

update_slider_text <- function(input, session, name_text, name_slider, min, max) {
  # Text was changed
  observeEvent(input[[name_text]],{
    input_N2 <- as.numeric(input[[name_text]])
    # Not a number, so change it to the slider
    if(is.na(input_N2)) {
      updateTextInput(
        session = session,
        inputId = name_text,
        value = input[[name_slider]]
      )
    }
    input_N2 <- pmin(max, pmax(min, input_N2))
    
    if(input_N2 != input[[name_slider]] &
       !is.na(input_N2) & input_N2 >= min & input_N2 <= max)
    {
      updateSliderInput(
        session = session,
        inputId = name_slider,
        value = input_N2
      )
    }
    
    # Also update text if pmin/pmax changed inputN2?
    if(input_N2 != input[[name_text]] & !is.na(input_N2)) {
      updateTextInput(
        session = session,
        inputId = name_text,
        value = input_N2
      ) 
    }
  })
  
  # Slider was changed
  observeEvent(input[[name_slider]],{
    if(as.numeric(input[[name_text]]) != input[[name_slider]] & !is.na(as.numeric(input[[name_text]])))
    {
      updateTextInput(
        session = session,
        inputId = name_text,
        value = input[[name_slider]]
      ) 
    }
  }, ignoreInit = T)
}

server <- function(input, output, session) {
  observeEvent(input$lowestexclude, {
    if (input$lowestexclude == "Lowest") {
      # withMathJax(updateSelectInput(session, "x_var", choices = c("r2", "Disease prevalence", "Number of embryos")))
      withMathJax(updateRadioButtons(session, "x_var", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos"), 
                                     choiceValues = c("r2", "Disease prevalence", "Number of embryos"), inline = T))
    }
    else {
      # withMathJax(updateSelectInput(session, "x_var", choices = c("r2", "Disease prevalence", "Number of embryos", "Quantile")))
      withMathJax(updateRadioButtons(session, "x_var", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos", "Quantile"),
                                     choiceValues = c("r2", "Disease prevalence", "Number of embryos", "Quantile"), inline = T))
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
    else if (input$x_var == "Quantile") {
      shinyjs::hide("q")
      shinyjs::show("N")
      shinyjs::show("K")
      shinyjs::show("r")
    }
  })
  
  # Should probably make it automatic?
  update_slider_text(input, session, "input_N2", "N2", 2, 1000)
  update_slider_text(input, session, "input_K2", "K2", 0.01, 1)
  update_slider_text(input, session, "input_r2", "r2", 0, 1)
  update_slider_text(input, session ,"input_q2", "q2", 0, 1)
  update_slider_text(input, session ,"input_qf2", "qf2", 0, 100)
  update_slider_text(input, session ,"input_qm2", "qm2", 0, 100)
  update_slider_text(input, session ,"input_h2", "h2", 0, 1)
  update_slider_text(input, session ,"input_samples", "samples", 5000, 300000)
  
  update_slider_text(input, session, "input_N_2", "N_2", 2, 100)
  update_slider_text(input, session, "input_rho", "rho", -0.99, 0.99)
  update_slider_text(input, session, "input_samples_2", "samples_2", 5000, 300000)
  update_slider_text(input, session ,"input_r2_1", "r2_1", 0, 1)
  update_slider_text(input, session ,"input_r2_2", "r2_2", 0, 1)
  update_slider_text(input, session ,"input_K_1", "K_1", 0, 100)
  update_slider_text(input, session ,"input_K_2", "K_2", 0, 1)

  output$distPlot <- renderPlot({
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
      x <- 1:20
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(input$r, input$K, n = x))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(r2 = input$r, K = input$K, q = input$q, n = x))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      
      # subtitle <- "Lowest strategy"
      x_lab <- "Number of embryos"
    }
    else if(selected_x == "Disease prevalence") {
      x <- seq(0.01, 1, length = 50)
      # y <- sapply(x, function(x) risk_reduction_lowest(input$r, x, n = input$N))
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(input$r, x, n = input$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * x
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(r2 = input$r, K = x, q = input$q, n = input$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * x
      }
      # subtitle <- "Lowest strategy"
      x_lab <- "Disease prevalence (log scale)"
    }
    else if(selected_x == "r2" || selected_x == "$$r^2$$") {
      x <- seq(0, 1, length = 50)
      # y <- sapply(x, function(x) risk_reduction_lowest(x, input$K, n = input$N))
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(x, input$K, n = input$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(x, input$K, input$q, n = input$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      # subtitle <- "Lowest strategy"
      # x_lab <- "PRS r^2"
      x_lab <- expression(r^2)
    }
    else if (selected_x == "Quantile") {
      # q
      x <- seq(0, 1, length = 50)
      y <- sapply(x, function(x) risk_reduction_exclude(input$r, input$K, x, n = input$N))
      # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      # subtitle <- "Exclude strategy"
      x_lab <- "Quantile to exclude"
    }
    else {
      print("Error!")
      print(selected_x)
    }
    
    ggplot(mapping = aes(x, y)) +
      geom_point() +
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
        plot.background = element_rect(fill = "#f5f5f5", color = NA)
      ) +
      scale_y_continuous(name = "Relative risk reduction",
                         sec.axis = sec_axis(~.*ifelse(selected_x == "K", x, input$K), name = "Absolute risk reduction"))
  })
  
  output$summary <- renderPrint({
    cat("<div class = \"alert alert-info\">")
    if (input$type2 == "Risk reduction") {
      temp <- risk_reduction_lowest(input$r2, K = input$K2, n = input$N2)
      if (input$lowestexclude2 != "Lowest") {
        temp <-
          risk_reduction_exclude(input$r2,
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
    }
    else if (input$type2 == "Conditional") {
      temp <-
        risk_reduction_lowest_conditional(
          input$r2,
          K = input$K2,
          n = input$N2,
          qf = input$qf2,
          qm = input$qm2,
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
            qf = input$qf2,
            qm = input$qm2,
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
    }
    else {
      if (input$r2 > input$h2) {
        updateNumericInput(session, "input_r2", value = input$h2)
        updateSliderInput(session, "r2", value = input$h2)
      }
      
      
      if (input$r2 <= input$h2) {
        if (input$lowestexclude2 != "Lowest") {
          temp <-
            risk_reduction_exclude_family_history2(
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
            risk_reduction_lowest_family_history2(
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
            "<p><u>Baseline risk</u>: <strong>%.4f(%.4f)</strong>\n</p>",
            temp[1],
            temp[5]
          )
        )
        cat(
          sprintf(
            "<p><u>Risk for specific strategy</u>: <strong>%.4f(%.4f)</strong>\n</p>",
            temp[2],
            temp[6]
          )
        )
        cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.4f</strong>\n</p>", temp[3]))
        cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.4f</strong>\n</p>", temp[4]))
        
        cat(sprintf("<p style = \"color:red\">Based on %d simulations, estimated standard deviation in parentheses.</p>", input$samples))
        if(temp[1] < temp[2]) {
          cat(sprintf("<p style = \"color:red\">Baseline risk is smaller than the strategy risk, so the sample size is probably too small.</p>", input$samples))
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
    cat(sprintf("<p><u>Relative risk reduction for disease 2</u>: <strong>%.4f</strong>\n</p>", temp[3]))
    cat(sprintf("<p><u>Absolute risk reduction for disease 2</u>: <strong>%.4f</strong>\n</p>", temp[4]))
    cat(sprintf("<p style = \"color:red\">Based on %d simulations.</p>", input$samples_2))
    cat("</div>")
  })
  
  output$about_text <- renderPrint({
    cat("This calcluator implements the methods from the paper Utility of polygenic embryo screening for
disease depends on the selection strategy.")
  })
}

# shinyApp(ui = ui, server = server)
# runApp(launch.browser = F)