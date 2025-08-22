###############################################################################
# app.R – Shiny demo, which simulates COVID-19 infections between 2020 and 2021
# Date  : 2025-07-25
###############################################################################

# Needed when setting-up shinylive, so packages are correctly included
if (FALSE) {
  library(shiny); library(shinyBS); library(ggplot2); library(data.table); 
  library(MicSim); library(munsell); library(colorspace); library(snowfall)
  library(rlecuyer); library(snow)
}

# loading libraries -------------------------------------------------------
library(shiny)
library(shinyBS)
library(MicSim)
library(data.table)
library(ggplot2)


# defining constants ------------------------------------------------------
startDate   <- 20200101   # yyyymmdd
endDate     <- 20211231   # yyyymmdd
simHorizon  <- c(startDate = startDate, endDate = endDate)
maxAge      <- 101
absStates   <- "dead"      # required by MicSim

# setting path to data
incMalePath   <- "./data/SurvStat_male_data.rds"
incFemalePath <- "./data/SurvStat_female_data.rds"


# function to run the sim -------------------------------------------------
# Defining a function to run the simulation,
# which includes the necessary pre and post processing
run_sim <- function(time_sick_days = 14,
                    effect_int     = 0,
                    prop_female    = 0.5,
                    date_int       = as.Date("2020-01-01"),
                    minage_start   = 18,
                    maxage_start   = 80,
                    N              = 1000) {
  
  # checking if the input is valid
  stopifnot(minage_start <= maxage_start,
            prop_female  >= 0, prop_female <= 1,
            effect_int   >= 0, effect_int <= 1,
            N > 0, time_sick_days > 0)
  
  # defining initial variables
  time_sick <- time_sick_days / 365
  N_female  <- trunc(N * (prop_female))
  N_male    <- N - N_female
  start_yr  <- as.Date(sprintf("%s-01-01", format(date_int, "%Y")))
  date_int  <- as.numeric(format(date_int, "%Y")) + as.numeric(date_int - start_yr) / 365.25
  
  # seed for reproducibility
  set.seed(9876)
  
  # Defining initial pop ----------------------------------------------------
  # setting reference data as jan 1st of the start year
  ref_date <- as.POSIXct(
    paste0(substr(startDate, 1, 4), "-01-01 00:00:00"),
    tz = "UTC"
  )
  
  # building function to return birthdays as strings in YYYYMMDD format
  init_bd  <- function(n) {
    ages <- runif(n, min = minage_start, max = maxage_start)
    format(as.Date(ref_date - ages * 365.25 * 24 * 3600), "%Y%m%d")
  }
  
  initPop <- data.frame(
    ID        = 1:N,
    birthDate = c(init_bd(N_female), init_bd(N_male)),
    initState = c(rep("Susceptible/Female", N_female),
                  rep("Susceptible/Male", N_male))
  )
  
  # loading the data --------------------------------------------------------
  inc_dat_m <- readRDS(incMalePath)
  inc_dat_f <- readRDS(incFemalePath)
  
  decYears  <- c(2020, inc_dat_m[[1]])
  inc_mat_m <- as.matrix(inc_dat_m[,-1])
  inc_mat_f <- as.matrix(inc_dat_f[,-1])
  ageBands  <- ncol(inc_mat_m)
  
  
  # defining transition rates -----------------------------------------------
  inc_rate_m <- function(age, calTime, duration) {
    col_index <- pmin(floor(age) + 1, ageBands)           # ages > 80 go to last band
    row_index <- findInterval(calTime, decYears)
    rate      <- ifelse(calTime >= date_int, 
                        (1 - effect_int) * inc_mat_m[cbind(row_index, col_index)],
                        inc_mat_m[cbind(row_index, col_index)])
    rate
  }
  assign("inc_rate_m", inc_rate_m, envir = .GlobalEnv)
  
  inc_rate_f <- function(age, calTime, duration) {
    col_index <- pmin(floor(age) + 1, ageBands)           # ages > 80 go to last band
    row_index <- findInterval(calTime, decYears)
    rate      <- ifelse(calTime >= date_int, 
                        (1 - effect_int) * inc_mat_f[cbind(row_index, col_index)],
                        inc_mat_f[cbind(row_index, col_index)])
    rate
  }
  assign("inc_rate_f", inc_rate_f, envir = .GlobalEnv)
  
  recovery_rate <- function(age, calTime, duration) ifelse(duration < time_sick, 0, Inf)
  assign("recovery_rate", recovery_rate, envir = .GlobalEnv)
  
  mortRates     <- function(age, calTime, duration) 0
  assign("mortRates", mortRates, envir = .GlobalEnv)
  
  
  # building transition matrix ----------------------------------------------
  health     <- c("Susceptible", "Infected", "Recovered")
  sex        <- c("Male", "Female")
  stateSpace <- expand.grid(health = health, sex = sex)
  
  TrMatrix_f <- cbind(c("Susceptible/Female->Infected/Female",
                        "Infected/Female->Recovered/Female"),
                      c("inc_rate_f", "recovery_rate"))
  TrMatrix_m <- cbind(c("Susceptible/Male->Infected/Male",
                        "Infected/Male->Recovered/Male"),
                      c("inc_rate_m", "recovery_rate"))
  allTransitions <- rbind(TrMatrix_f, TrMatrix_m)
  absTransitions <- cbind("dead", "mortRates")
  
  transitionMatrix <- buildTransitionMatrix(allTransitions = allTransitions,
                                            absTransitions = absTransitions,
                                            stateSpace      = stateSpace)
  
  
  # running the simulation --------------------------------------------------
  # ensuring these are all available in global environment
  assign("stateSpace",     stateSpace,     envir = .GlobalEnv)
  assign("allTransitions", allTransitions, envir = .GlobalEnv)
  assign("absTransitions", absTransitions, envir = .GlobalEnv)
  assign("simHorizon",     simHorizon,     envir = .GlobalEnv)
  
  pop <- micSim(initPop          = initPop,
                transitionMatrix = transitionMatrix,
                absStates        = absStates,
                maxAge           = maxAge,
                simHorizon       = simHorizon)
  
  # handling case where no transitions
  if (nrow(pop) == N & sum(is.na(pop$From)) == N) {
    # building date range
    start_d <- as.IDate(as.character(simHorizon["startDate"]), "%Y%m%d")
    end_d   <- as.IDate(as.character(simHorizon["endDate"]),   "%Y%m%d")
    dates   <- seq(start_d, end_d, by = 1L)
    n_days  <- length(dates)
    
    # assigning metrics by sex
    metrics <- data.table(
      sex                = rep(c("Female","Male", "Overall"), each = n_days),
      date               = rep(dates, times = 3L),
      active_inf_14d     = 0L,
      cum_recovered      = 0L,
      susceptible        = c(rep.int(N_female, n_days), 
                             rep.int(N_male, n_days),
                             rep.int(N, n_days)),
      N                  = c(rep.int(N_female, n_days), 
                             rep.int(N_male, n_days),
                             rep.int(N, n_days)),
      active_inf_14d_pct = 0L,
      cum_recovered_pct  = 0L,
      susceptible_pct    = 100L
    )
    
    return(metrics)
  }
  
  
  # processing the output ---------------------------------------------------
  long <- as.data.table(convertToLongFormat(pop))
  long[, `:=`(start = as.IDate(Tstart, "%Y%m%d"),
              stop  = as.IDate(Tstop,  "%Y%m%d"))]
  
  all_dates <- CJ(sex = unique(long$sex), date = seq(min(long$start), max(long$stop), by = 1))
  counts    <- long[, .(count = uniqueN(ID)), by = .(sex, start, health)]
  
  # getting number of infections within 14 day window
  inf <- merge(all_dates,
               counts[health == "Infected", .(sex, date = start, new_inf = count)],
               by = c("sex", "date"), all.x = TRUE)[
                 , new_inf := fifelse(is.na(new_inf), 0L, new_inf)][
                   , active_inf_14d := frollsum(new_inf, 14, align = "right", fill = 0),
                   by = sex]
  
  # getting number recovered
  rec <- merge(all_dates,
               counts[health == "Recovered", .(sex, date = start, new_rec = count)],
               by = c("sex", "date"), all.x = TRUE)[
                 , new_rec := fifelse(is.na(new_rec), 0L, new_rec)][
                   , cum_recovered := cumsum(new_rec), by = sex]
  
  # getting number susceptible
  init_sus <- counts[start == min(long$start) & health == "Susceptible",
                     .(sex, init_sus = count)]
  sus <- merge(inf[, .(sex, date, new_inf)], init_sus, by = "sex")[
    , susceptible := init_sus - cumsum(new_inf), by = sex]
  
  # combining all the metrics
  metrics_sex <- Reduce(function(x, y) merge(x, y, by = c("sex", "date")),
                        list(inf[, .(sex, date, active_inf_14d)],
                             rec[, .(sex, date, cum_recovered)],
                             sus[, .(sex, date, susceptible)]))
  
  metric_cols <- c("active_inf_14d", "cum_recovered", "susceptible")
  N_sex       <- data.table(sex = c("Female", "Male", "Overall"), N = c(N_female, N_male, N))
  
  overall <- metrics_sex[, lapply(.SD, sum), by = date, .SDcols = metric_cols][
    , sex := "Overall"]
  
  metrics <- rbind(metrics_sex, overall)
  metrics <- merge(metrics, N_sex, by = "sex")
  
  metrics[, `:=`(active_inf_14d_pct = active_inf_14d / N * 100,
                 cum_recovered_pct  = cum_recovered  / N * 100,
                 susceptible_pct    = susceptible    / N * 100)]
  metrics
}


# Shiny UI ----------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .container-fluid { max-width: 1200px; }
    .btn { border-radius: 8px; }
    .form-control, .irs { border-radius: 8px; }
    .shiny-output-error-validation { color: #b00020; }
  "))),
  titlePanel("Simulation of COVID-19 Infections in Germany (2020–2021)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      numericInput("N", "Sample size (N)",
                   value = 1000, min = 1000, max = 10000, step = 500),
      # sliderInput ("prop_female", "Proportion female (%)",
      # min = 0,  max = 100, value = 50, step = 1),
      
      sliderInput("minage_start", "Minimum age at start", 
                  min = 0,  max = 100, value = 18),
      sliderInput("maxage_start", "Maximum age at start", 
                  min = 0,  max = 100, value = 80),
      
      numericInput("time_sick", "Duration of illness (days)",
                   value = 14,  min = 1,  max = 60),
      
      sliderInput("effect_int", "Intervention effectiveness (%)",
                  min = 0,  max = 100, value = 0, step = 1),
      bsPopover(
        id        = "effect_int",
        title     = "What does this do?",
        content   = paste(
          "Reduces all incidence rates by this percent from the chosen start date.",
          "For example, 25% means the incidence rates are all multiplied by 0.75."
        ),
        placement = "right",
        trigger   = "hover"
      ),
      
      conditionalPanel(
        condition = "input.effect_int > 0",
        dateInput(inputId = "date_int", 
                  label = "Intervention starts", 
                  value = "2020-01-01", min = "2020-01-01", max = "2021-12-31", 
                  format = "yyyy-mm-dd")
      ),
      
      div(
        class = "form-group shiny-input-container mb-3",  # matches Shiny input wrappers (BS3/BS5+bslib)
        tags$label("Display options", class = "control-label form-label"),
        checkboxInput("by_sex", "Break down by sex", FALSE),
        checkboxGroupInput("series", "Show series",
                           choices  = c("Infected (14-day active)" = "inf",
                                        "Recovered (cumulative)"   = "rec",
                                        "Susceptible"              = "sus"),
                           selected = c("inf","rec")
        )
      ),
      #checkboxInput("by_sex", "Break down by sex", FALSE),
      
      # checkboxGroupInput("series", "Show series", 
      #                    choices = c("Infected (14‑day active)" = "inf", 
      #                                "Recovered (cumulative)"   = "rec",
      #                                "Susceptible"              = "sus"),
      #                    selected = c("inf","rec")
      # ),
      
      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    
    mainPanel(
      width = 9,
      fluidRow(
        column(4, uiOutput("kpi_peak")),
        column(4, uiOutput("kpi_peak_date")),
        column(4, uiOutput("kpi_recovered"))
      ),
      plotOutput("epiPlot", height = "600px"),
      downloadButton("dl_csv", "Download metrics (CSV)", class = "btn btn-outline-secondary me-2"),
      downloadButton("dl_png", "Download plot (PNG)", class = "btn btn-outline-secondary")
    )
  )
)


# Shiny server ------------------------------------------------------------
server <- function(input, output, session) {
  
  # keep age sliders consistent
  observeEvent(input$minage_start, {
    if (input$minage_start > input$maxage_start)
      updateSliderInput(session, "maxage_start", value = input$minage_start)
  })
  observeEvent(input$maxage_start, {
    if (input$maxage_start < input$minage_start)
      updateSliderInput(session, "minage_start", value = input$maxage_start)
  })
  
  # storing results of last simulation run
  metricsData <- eventReactive(input$run, {
    
    # input validation
    validate(
      need(input$minage_start <= input$maxage_start, "Minimum age must be ≤ maximum age."),
      need(input$N > 0,                          "Sample size must be positive."),
      need(input$time_sick > 0,                  "Duration of illness must be positive."),
      need(file.exists(incMalePath) && file.exists(incFemalePath),
           "Incidence data not found in the data folder.")
    )
    
    withProgress(message = "Running simulation...", value = 0.1, {
      # call helper
      metrics <- run_sim(time_sick_days = input$time_sick,
                         effect_int     = input$effect_int / 100,
                         #prop_female    = 0.5, #input$prop_female / 100,
                         minage_start   = input$minage_start,
                         maxage_start   = input$maxage_start,
                         N              = input$N,
                         date_int       = input$date_int)
      incProgress(0.9)
      metrics
    })
  })
  
  kpi_box <- function(title, value) {
    div(class="p-3 border rounded mb-3",
        tags$div(class="text-muted small", title),
        tags$div(class="h4 mb-0", value))
  }
  
  kpi_dual <- function(title, female_value, male_value) {
    div(class="p-3 border rounded mb-3",
        tags$div(class="text-muted small", title),
        div(class="d-flex justify-content-between gap-3",
            div(class="flex-fill",
                tags$div(class="small text-muted", "Female"),
                tags$div(class="h4 mb-0", female_value)
            ),
            div(class="flex-fill",
                tags$div(class="small text-muted", "Male"),
                tags$div(class="h4 mb-0", male_value)
            )
        )
    )
  }
  
  output$kpi_peak <- renderUI({
    dat <- metricsData()
    req(nrow(dat))
    
    if (isTRUE(input$by_sex)) {
      df <- dat[sex == "Female"]; dm <- dat[sex == "Male"]
      kpi_dual(
        "Peak infected (14-day active)",
        sprintf("%0.1f%%", max(df$active_inf_14d_pct, na.rm = TRUE)),
        sprintf("%0.1f%%", max(dm$active_inf_14d_pct, na.rm = TRUE))
      )
    } else {
      d <- dat[sex == "Overall"]
      kpi_box("Peak infected (14-day active)",
              sprintf("%0.1f%%", max(d$active_inf_14d_pct, na.rm = TRUE)))
    }
  })
  
  output$kpi_peak_date <- renderUI({
    dat <- metricsData()
    req(nrow(dat))
    
    if (isTRUE(input$by_sex)) {
      df <- dat[sex == "Female"]; dm <- dat[sex == "Male"]
      idx_f <- which.max(df$active_inf_14d_pct)
      idx_m <- which.max(dm$active_inf_14d_pct)
      kpi_dual(
        "Peak date",
        format(df$date[idx_f], "%b %d, %Y"),
        format(dm$date[idx_m], "%b %d, %Y")
      )
    } else {
      d <- dat[sex == "Overall"]; idx <- which.max(d$active_inf_14d_pct)
      kpi_box("Peak date", format(d$date[idx], "%b %d, %Y"))
    }
  })
  
  output$kpi_recovered <- renderUI({
    dat <- metricsData()
    req(nrow(dat))
    
    if (isTRUE(input$by_sex)) {
      df <- dat[sex == "Female"]; dm <- dat[sex == "Male"]
      kpi_dual(
        "Recovered by end",
        sprintf("%0.1f%%", tail(df$cum_recovered_pct, 1)),
        sprintf("%0.1f%%", tail(dm$cum_recovered_pct, 1))
      )
    } else {
      d <- dat[sex == "Overall"]
      kpi_box("Recovered by end", sprintf("%0.1f%%", tail(d$cum_recovered_pct, 1)))
    }
  })
  
  plotObj <- reactive({
    dat <- metricsData()
    req(nrow(dat))
    
    plot_dat <- if (isTRUE(input$by_sex)) dat[sex != "Overall"] else dat[sex == "Overall"]
    
    y_map <- list(
      inf = "active_inf_14d_pct",
      rec = "cum_recovered_pct",
      sus = "susceptible_pct"
    )
    chosen <- unlist(y_map[input$series], use.names = FALSE)
    req(length(chosen) > 0)
    
    long <- data.table::melt(
      plot_dat,
      id.vars = c("sex","date"),
      measure.vars = chosen,
      variable.name = "metric",
      value.name = "pct"
    )
    metric_labels <- c(active_inf_14d_pct = "Infected (14‑day active)",
                       cum_recovered_pct  = "Recovered (cumulative)",
                       susceptible_pct    = "Susceptible")
    long[, metric := factor(metric, levels = names(metric_labels), labels = unname(metric_labels))]
    
    p <- ggplot(long, aes(x = date, y = pct, colour = metric)) +
      geom_line(linewidth = 0.9) +
      { if (isTRUE(input$by_sex)) facet_wrap(~sex, ncol = 1) else NULL } +
      scale_y_continuous(labels = scales::label_percent(accuracy = 0.1, scale = 1)) +
      labs(x = "Date", y = "Percent of cohort", colour = NULL) +
      theme_light(base_size = 14) +
      theme(legend.position = "top", panel.grid.minor = element_blank())
    
    if (!is.null(input$effect_int) && input$effect_int > 0 && !is.null(input$date_int)) {
      p <- p +
        geom_vline(xintercept = as.Date(input$date_int), linetype = 2) +
        annotate("label", x = as.Date(input$date_int), y = Inf, vjust = 1,
                 label = "Intervention", size = 3)
    }
    p
  })
  
  # render plot
  output$epiPlot <- renderPlot({ plotObj() })
  
  output$dl_csv <- downloadHandler(
    filename = function() sprintf("simulation_metrics_%s.csv", Sys.Date()),
    content = function(file) data.table::fwrite(metricsData(), file)
  )
  output$dl_png <- downloadHandler(
    filename = function() sprintf("simulation_plot_%s.png", Sys.Date()),
    content = function(file) {
      ggplot2::ggsave(filename = file, plot = plotObj(), width = 10, height = 6, dpi = 150)
    }
  )

}

# launch ------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)
app
