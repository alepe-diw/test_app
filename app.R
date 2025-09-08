###############################################################################
# app.R – Shiny demo, which simulates COVID-19 infections between 2020 and 2021
# Date  : 2025-09-08
###############################################################################

# needed when setting-up shinylive, so packages are correctly included
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
# defining a function to run the simulation,
# which includes the necessary pre and post processing
run_sim <- function(time_sick_days = 14,
                    effect_int     = 0,
                    prop_female    = .5,
                    date_int       = as.Date("2020-01-01"),
                    minage_start   = 18,
                    maxage_start   = 80,
                    N              = 1000) {
  
  # checking if the input is valid
  stopifnot(prop_female  >= 0, prop_female <= 1,
            effect_int   >= 0, effect_int  <= 1,
            N >= 100, time_sick_days > 0)
  
  # defining initial variables
  time_sick <- time_sick_days / 365
  N_female  <- trunc(N * (prop_female))
  N_male    <- N - N_female
  start_yr  <- as.Date(sprintf("%s-01-01", format(date_int, "%Y")))
  date_int  <- as.numeric(format(date_int, "%Y")) + 
    as.numeric(date_int - start_yr) / 365.25
  
  # seed for reproducibility
  set.seed(9876)
  
  
  # defining initial pop ----------------------------------------------------
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
    
    rate <- ifelse(calTime >= date_int,
                   (1 - effect_int) * inc_mat_m[cbind(row_index, col_index)],
                   inc_mat_m[cbind(row_index, col_index)])
    rate
  }
  assign("inc_rate_m", inc_rate_m, envir = .GlobalEnv)
  
  inc_rate_f <- function(age, calTime, duration) {
    col_index <- pmin(floor(age) + 1, ageBands)           # ages > 80 go to last band
    row_index <- findInterval(calTime, decYears)
    
    rate <- ifelse(calTime >= date_int,
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
      active_inf         = 0L,
      active_inf_7d      = 0L,
      cum_recovered      = 0L,
      susceptible        = c(rep.int(N_female, n_days), 
                             rep.int(N_male, n_days),
                             rep.int(N, n_days)),
      N                  = c(rep.int(N_female, n_days), 
                             rep.int(N_male, n_days),
                             rep.int(N, n_days)),
      active_inf_pct     = 0L,
      active_inf_7d_pct  = 0L,
      cum_recovered_pct  = 0L,
      susceptible_pct    = 100L
    )
    
    return(metrics)
  }
  
  
  # processing the output ---------------------------------------------------
  long <- as.data.table(convertToLongFormat(pop))
  long[, `:=`(start = fix_yyyymmdd(Tstart),
              stop  = fix_yyyymmdd(Tstop))]
  
  all_dates <- CJ(sex = unique(long$sex), date = seq(min(long$start), max(long$stop), by = 1))
  counts    <- long[, .(count = uniqueN(ID)), by = .(sex, start, health)]
  
  # getting number of infections within a window determined by time_sick_days
  inf <- merge(all_dates,
               counts[health == "Infected", .(sex, date = start, new_inf = count)],
               by = c("sex", "date"), all.x = TRUE)[
                 , new_inf := fifelse(is.na(new_inf), 0L, new_inf)][
                   , `:=`(
                     active_inf    = frollsum(new_inf, time_sick_days, align = "right", fill = 0),
                     active_inf_7d = frollsum(new_inf, 7, align = "right", fill = 0)
                   ),
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
                        list(inf[, .(sex, date, active_inf, active_inf_7d)],
                             rec[, .(sex, date, cum_recovered)],
                             sus[, .(sex, date, susceptible)]))
  
  metric_cols <- c("active_inf", "active_inf_7d", "cum_recovered", "susceptible")
  N_sex       <- data.table(sex = c("Female", "Male", "Overall"), N = c(N_female, N_male, N))
  
  overall <- metrics_sex[, lapply(.SD, sum), by = date, .SDcols = metric_cols][
    , sex := "Overall"]
  
  metrics <- rbind(metrics_sex, overall)
  metrics <- merge(metrics, N_sex, by = "sex")
  
  metrics[, `:=`(active_inf_pct    = active_inf    / N * 100,
                 active_inf_7d_pct = active_inf_7d / N * 100,
                 cum_recovered_pct = cum_recovered / N * 100,
                 susceptible_pct   = susceptible   / N * 100)]
  metrics
}


# helper functions --------------------------------------------------------
# outputs key indicators from the sim
summary_indicators <- function(dt) {
  list(
    peak_inf_pct  = if (nrow(dt)) max(dt$active_inf_7d_pct, na.rm = TRUE) else NA_real_,
    peak_inf_date = if (nrow(dt)) dt$date[which.max(dt$active_inf_7d_pct)] else as.Date(NA),
    end_rec_pct   = if (nrow(dt)) tail(dt$cum_recovered_pct, 1) else NA_real_,
    end_sus_pct   = if (nrow(dt)) tail(dt$susceptible_pct, 1) else NA_real_
  )
}

# formatting percentages
fmt_pct <- function(x, d=1) ifelse(is.finite(x), sprintf(paste0("%0.", d, "f%%"), x), "—")

# assigning colors and line types to different indicators. 
# using SOEP colors for now
metric_colors <- c(
  active_inf_pct    = "#5c5c5c",
  active_inf_7d_pct = "#ae393f",
  cum_recovered_pct = "#00786b",
  susceptible_pct   = "#9467bd" 
)
linetype_vals <- c("Baseline" = "solid", "Intervention" = "longdash")

# ui helpers for making nicer boxes
indicator_box <- function(title, value) {
  div(class = "p-3 border rounded mb-3",
      tags$div(class = "text-muted small", title),
      tags$div(class = "h4 mb-0", value))
}
indicator_dual <- function(title, female_value, male_value) {
  div(class = "p-3 border rounded mb-3",
      tags$div(class = "text-muted small", title),
      div(class = "d-flex justify-content-between gap-3",
          div(class = "flex-fill",
              tags$div(class = "small text-muted", "Female"),
              tags$div(class = "h4 mb-0", female_value)),
          div(class = "flex-fill",
              tags$div(class = "small text-muted", "Male"),
              tags$div(class = "h4 mb-0", male_value))))
}

# drops overall group when the user wants to see results by sex
filter_sex <- function(dt, by_sex) {
  if (isTRUE(by_sex)) dt[sex != "Overall"] else dt[sex == "Overall"]
}

# fixes issues where dates in yyyymmdd format are invalid, e.g., 20210229
# any day past the end of the month is assumed to have occured on the last day
# of the month
fix_yyyymmdd <- function(x) {
  # trying the simple conversion first and returns if no errors
  out <- data.table::as.IDate(x, "%Y%m%d")
  bad <- is.na(out) & !is.na(x)
  if (!any(bad)) return(out)
  
  # if there are errors, fixing them
  # this first part mimics the existing code in getInDateFormat
  xi <- as.integer(x[bad])
  y  <- xi %/% 10000L
  m  <- (xi %/% 100L) %% 100L
  d  <- xi %% 100L
  
  # safeguarding against month being 13. If the month is 13, interprets
  # it as January of the next year
  y <- y + (m == 13L)
  m[m == 13L] <- 1L
  
  # handling incorrect number of days in a month
  # determining max number of days in a month (with leap years)
  dim_month <- c(31L,28L,31L,30L,31L,30L,31L,31L,30L,31L,30L,31L)
  maxd      <- dim_month[m]
  is_leap   <- (y %% 400L == 0L) | ((y %% 4L == 0L) & (y %% 100L != 0L))
  feb_idx   <- (m == 2L)
  if (any(feb_idx)) maxd[feb_idx] <- 28L + as.integer(is_leap[feb_idx])
  
  # ensuring the day value is between 1 and max number of day
  d <- pmax.int(1L, pmin.int(d, maxd))
  
  # returning a date
  out[bad] <- data.table::as.IDate(sprintf("%04d%02d%02d", y, m, d), "%Y%m%d")
  out
}

# outputs the values from summary_indicators() to the ui
output_indicators_row <- function(df, label = "", title_prefix = "", days = 14) {
  ind <- summary_indicators(df)
  suffix <- if (nzchar(label)) paste0(" — ", label) else ""
  fluidRow(
    column(3, indicator_box(
      paste0(title_prefix, "Peak infected (7-day total)", suffix),
      #sprintf("%sPeak infected (%d-day active)%s", title_prefix, days, suffix),
      fmt_pct(ind$peak_inf_pct)
    )),
    column(3, indicator_box(
      paste0(title_prefix, "Peak date", suffix),
      format(ind$peak_inf_date, "%b %d, %Y")
    )),
    column(3, indicator_box(
      paste0(title_prefix, "Recovered by end", suffix),
      fmt_pct(ind$end_rec_pct)
    )),
    column(3, indicator_box(
      paste0(title_prefix, "Susceptible by end", suffix),
      fmt_pct(ind$end_sus_pct)
    ))
  )
}


# shiny UI ----------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Simulation of COVID-19 Infections in Germany (2020–2021)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Cohort & disease"),
      numericInput("N", "Sample size (N)",
                   value = 1000, min = 1000,
                   max = 30000, step = 1000),
      sliderInput("minage_start", "Minimum age at start",
                  min = 0,  max = 100, value = 18),
      sliderInput("maxage_start", "Maximum age at start",
                  min = 0,  max = 100, value = 80),
      numericInput("time_sick", "Duration of illness (days)",
                   value = 14,  min = 1,  max = 60, step = 1),
      
      tags$hr(),
      h4("Intervention scenario"),
      sliderInput("effect_int", "Intervention effectiveness (%)",
                  min = 0,  max = 99, value = 0, step = 1),
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
      conditionalPanel(
        condition = "input.effect_int > 0",
        dateInput("date_int", "Intervention starts", value = "2021-01-01",
                  min = "2020-01-01", max = "2021-12-31", format = "yyyy-mm-dd"),
        checkboxInput("compare_baseline", "Compare intervention against baseline", FALSE)
      ),
      
      tags$hr(),
      h4("Display options"),
      checkboxInput("by_sex", "Break down by sex", FALSE),
      checkboxGroupInput("series", "Show series",
                         choices = c(
                           "Infected (7-day total)"               = "inf7",
                           "Infected (active, user-set duration)" = "inf",
                           "Recovered (cumulative)"               = "rec",
                           "Susceptible"                          = "sus"
                         ),
                         selected = c("inf7", "rec")),
      
      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      uiOutput("indicators_summary"),
      conditionalPanel(
        condition = "input.effect_int > 0 && input.compare_baseline == true",
        uiOutput("indicator_diff")
      ),
      plotOutput("epiPlot", height = "600px"),
      div(class = "mt-3",
          downloadButton("dl_csv", "Download metrics (CSV)", class = "btn btn-outline-secondary me-2"),
          downloadButton("dl_png", "Download plot (PNG)", class = "btn btn-outline-secondary"))
    )
  )
)


# Shiny server ------------------------------------------------------------
server <- function(input, output, session) {
  
  
  # validation --------------------------------------------------------------
  # keep age sliders consistent
  observeEvent(input$minage_start, {
    if (input$minage_start > input$maxage_start)
      updateSliderInput(session, "maxage_start", value = input$minage_start)
  })
  observeEvent(input$maxage_start, {
    if (input$maxage_start < input$minage_start)
      updateSliderInput(session, "minage_start", value = input$maxage_start)
  })
  
  
  # running the simulation(s) -----------------------------------------------
  sims <- eventReactive(input$run, {
    validate(
      need(input$N >= 1000, "Sample size must be at least 1,000 (larger values are recommended)."),
      need(input$time_sick > 0, "Duration of illness must be positive."),
      need((input$time_sick %% 1) == 0, "Duration of illness must be a whole number."),
      need(file.exists(incMalePath) && file.exists(incFemalePath),
           "Incidence data not found in the data folder.")
    )
    
    compare_mode <- (input$effect_int > 0) && isTRUE(input$compare_baseline)
    
    # if comparing baseline and intervention, run both
    if (compare_mode) {
      withProgress(message = "Running baseline and intervention...", value = 0.1, {
        # run baseline scenario 
        base <- run_sim(
          time_sick_days = input$time_sick,
          effect_int     = 0,
          minage_start   = input$minage_start,
          maxage_start   = input$maxage_start,
          N              = input$N
        )
        incProgress(0.5, detail = "Baseline done…")
        
        # run intervention scenario 
        interv <- run_sim(
          time_sick_days = input$time_sick,
          effect_int     = input$effect_int / 100,
          minage_start   = input$minage_start,
          maxage_start   = input$maxage_start,
          N              = input$N,
          date_int       = input$date_int
        )
        incProgress(0.9, detail = "Intervention done…")
        list(mode = "compare", baseline = base, intervention = interv, time_sick = input$time_sick)
      })
    } 
    else {
      withProgress(message = "Running scenario…", value = 0.3, {
        if (input$effect_int > 0 && !is.null(input$date_int)) {
          scen <- run_sim(
            time_sick_days = input$time_sick,
            effect_int     = input$effect_int / 100,
            minage_start   = input$minage_start,
            maxage_start   = input$maxage_start,
            N              = input$N,
            date_int       = input$date_int
          )
        } else {
          scen <- run_sim(
            time_sick_days = input$time_sick,
            effect_int     = input$effect_int / 100,
            minage_start   = input$minage_start,
            maxage_start   = input$maxage_start,
            N              = input$N
          )
        }
        incProgress(0.9)
        list(mode = "single", scenario = scen, time_sick = input$time_sick)
      })
    }
  })
  
  
  # summary figures ---------------------------------------------------------
  # getting key indicators
  output$indicators_summary <- renderUI({
    s <- sims(); req(s)
    days <- s$time_sick
    
    if (s$mode == "single") {
      dat <- filter_sex(s$scenario, input$by_sex)
      prefix <- ""                    # no "Scenario:" text
    } else {
      dat <- filter_sex(s$baseline, input$by_sex)
      prefix <- "Baseline: "
    }
    
    if (isTRUE(input$by_sex)) {
      df <- dat[sex == "Female"]; dm <- dat[sex == "Male"]
      tagList(
        output_indicators_row(df, "Female", prefix, days = days),
        output_indicators_row(dm, "Male",   prefix, days = days)
      )
    } else {
      output_indicators_row(dat, "", prefix, days = days)
    }
  })
  
  # calculating differences when compare is on
  output$indicator_diff <- renderUI({
    s <- sims(); req(s$mode == "compare", s$baseline, s$intervention)
    
    make_boxes <- function(slice_label, db, di) {
      kb <- summary_indicators(db); ki <- summary_indicators(di)
      d_peak     <- ki$peak_inf_pct - kb$peak_inf_pct
      d_end_rec  <- ki$end_rec_pct  - kb$end_rec_pct
      d_end_sus  <- ki$end_sus_pct  - kb$end_sus_pct
      fluidRow(
        column(4, indicator_box(paste0("Δ Peak infected (pp) — ", slice_label),
                                sprintf("%+0.1f pp (%s vs %s)", d_peak, fmt_pct(ki$peak_inf_pct), fmt_pct(kb$peak_inf_pct)))),
        column(4, indicator_box(paste0("Δ Recovered by end (pp) — ", slice_label),
                                sprintf("%+0.1f pp (%s vs %s)", d_end_rec,  fmt_pct(ki$end_rec_pct),  fmt_pct(kb$end_rec_pct)))),
        column(4, indicator_box(paste0("Δ Susceptible by end (pp) — ", slice_label),
                                sprintf("%+0.1f pp (%s vs %s)", d_end_sus,  fmt_pct(ki$end_sus_pct),  fmt_pct(kb$end_sus_pct))))
      )
    }
    
    if (isTRUE(input$by_sex)) {
      dbf <- s$baseline[sex=="Female"]; dimf <- s$intervention[sex=="Female"]
      dbm <- s$baseline[sex=="Male"];   dimm <- s$intervention[sex=="Male"]
      tagList(
        h5("Intervention vs baseline: indicator changes (percentage points)"),
        make_boxes("Female", dbf, dimf),
        make_boxes("Male",   dbm, dimm)
      )
    } else {
      db <- s$baseline[sex=="Overall"]; di <- s$intervention[sex=="Overall"]
      tagList(
        h5("Intervention vs baseline: indicator changes (percentage points)"),
        make_boxes("Overall", db, di)
      )
    }
  })
  
  
  # plots -------------------------------------------------------------------
  # assigning labels to indicators
  metric_labels <- reactive({
    days <- req(sims())$time_sick
    c(
      active_inf_pct    = sprintf("Infected (%d-day active)", days),
      active_inf_7d_pct = "Infected (7-day total)",
      cum_recovered_pct = "Recovered (cumulative)",
      susceptible_pct   = "Susceptible"
    )
  })
  
  # generating the plots
  plotObj <- reactive({
    s <- sims(); req(s)
    
    y_map <- list(inf  = "active_inf_pct",
                  inf7 = "active_inf_7d_pct",
                  rec  = "cum_recovered_pct",
                  sus  = "susceptible_pct")
    chosen <- unlist(y_map[input$series], use.names = FALSE); req(length(chosen) > 0)
    
    # getting labels for plots
    labs_map <- metric_labels()
    
    # extracting simulation output
    if (s$mode == "compare") {
      base <- if (isTRUE(input$by_sex)) s$baseline[sex != "Overall"] else s$baseline[sex == "Overall"]
      intv <- if (isTRUE(input$by_sex)) s$intervention[sex != "Overall"] else s$intervention[sex == "Overall"]
      base[, scenario := "Baseline"]; intv[, scenario := "Intervention"]
      wide <- rbind(base, intv, use.names = TRUE)
    } else {
      wide <- if (isTRUE(input$by_sex)) s$scenario[sex != "Overall"] else s$scenario[sex == "Overall"]
    }
    
    # pivoting long for ggplot
    long <- data.table::melt(
      wide,
      id.vars      = intersect(names(wide), c("sex", "date", "scenario")),
      measure.vars = chosen,
      variable.name = "metric",
      value.name    = "pct"
    )
    long[, metric := factor(metric, levels = names(labs_map))]
    
    p <- ggplot(long, aes(x = date, y = pct, colour = metric)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(
        values = metric_colors[names(labs_map)],
        breaks = names(labs_map),
        labels = unname(labs_map),
        name = "Metric"
      ) +
      scale_y_continuous(labels = scales::comma_format(accuracy = 0.1)) +
      labs(x = "Date", y = "Percentage of total cohort", colour = NULL) +
      theme_light(base_size = 14) +
      theme(legend.position = "top")
    
    # conditional things added to figure based on the users selections
    if (s$mode == "compare") {
      p <- p +
        aes(linetype = scenario) +
        guides(
          colour   = guide_legend(order = 1, override.aes = list(linetype = "solid", linewidth = 1.2)),
          linetype = guide_legend(title = "Scenario", order = 2)
        )
    }
    
    if (isTRUE(input$by_sex)) {
      p <- p + facet_wrap(~sex, ncol = 1)
    }
    
    if ((input$effect_int > 0) && !is.null(input$date_int)) {
      p <- p + geom_vline(xintercept = as.Date(input$date_int), linetype = 2)
    }
    
    p
  })
  
  output$epiPlot <- renderPlot({ plotObj() })
  
  
  # download output ---------------------------------------------------------
  output$dl_csv <- downloadHandler(
    filename = function() sprintf("simulation_metrics_%s.csv", Sys.Date()),
    content = function(file) {
      s <- sims(); req(s)
      if (s$mode=="compare") {
        b <- data.table::copy(s$baseline)[, scenario := "Baseline"]
        i <- data.table::copy(s$intervention)[, scenario := "Intervention"]
        out <- rbind(b, i, use.names = TRUE)
      } else {
        out <- data.table::copy(s$scenario)[, scenario := if (input$effect_int > 0) "Intervention" else "Baseline"]
      }
      data.table::fwrite(out, file)
    }
  )
  output$dl_png <- downloadHandler(
    filename = function() sprintf("simulation_plot_%s.png", Sys.Date()),
    content = function(file) {
      p <- plotObj()
      ggplot2::ggsave(filename = file, plot = p, width = 10, height = 6, dpi = 150)
    }
  )
}


# launch ------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)
app