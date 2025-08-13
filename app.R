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
                    prop_female    = .5,
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
    (1 - effect_int) * inc_mat_m[cbind(row_index, col_index)]
  }
  assign("inc_rate_m", inc_rate_m, envir = .GlobalEnv)
  
  inc_rate_f <- function(age, calTime, duration) {
    col_index <- pmin(floor(age) + 1, ageBands)           # ages > 80 go to last band
    row_index <- findInterval(calTime, decYears)
    (1 - effect_int) * inc_mat_f[cbind(row_index, col_index)]
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
  titlePanel("Simulation of COVID-19 Infections in Germany (2020–2021)"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 numericInput("N", "Sample size (N)",
                              value = 1000, min = 1000,
                              max = 10000, step = 500),
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
                     "This applies a percentage reduction to all the incidence rates.",
                     "For example, 25% means the incidence rates are all multiplied by 0.75."
                   ),
                   placement = "right",
                   trigger   = "hover"
                 ),
                 actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(width = 9,
              plotOutput("epiPlot", height = "600px")
    )
  )
)


# Shiny server ------------------------------------------------------------
server <- function(input, output, session) {
  
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
                         N              = input$N)
      incProgress(0.9)
      metrics
    })
  })
  
  # plot output
  output$epiPlot <- renderPlot({
    metrics_dat <- metricsData()
    req(metrics_dat)
    
    plot_dat <- metrics_dat[sex == "Overall", ]
    
    ggplot(plot_dat, aes(x = date)) +
      geom_line(aes(y = active_inf_14d_pct, colour = "Infected")) +
      geom_line(aes(y = cum_recovered_pct,  colour = "Recovered")) +
      #facet_wrap(.~sex) +
      scale_y_continuous(labels = scales::comma_format(accuracy = 0.1)) +
      labs(x = "Date", y = "Percentage of total cohort", colour = NULL) +
      theme_minimal(base_size = 14)
  })
}

# launch ------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)
app