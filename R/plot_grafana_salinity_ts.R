GRAFANA_BASE <- "https://mbon-dashboards.marine.usf.edu"
GRAFANA_DASHBOARD_UID <- "cfraol6xev5kwd"
GRAFANA_DASHBOARD_SLUG <- "south-florida-salinity"
GRAFANA_SALINITY_PANEL_ID <- 33L
GRAFANA_INFLUX_DATASOURCE_UID <- "P3C6603E967DC8568"
GRAFANA_DEFAULT_FROM <- "1720965498428"
GRAFANA_DEFAULT_TO <- "1784037498428"

grafana_buoy_label <- function(station_info) {
  short <- sub("f1$", "", station_info$name)
  sprintf("%s (%s)", station_info$long_name, short)
}

grafana_regex_buoy <- function(buoy_label) {
  gsub("([][.+*?^$(){}|\\\\])", "\\\\\\1", buoy_label, perl = TRUE)
}

grafana_solo_panel_url <- function(
    buoy_label,
    from = GRAFANA_DEFAULT_FROM,
    to = GRAFANA_DEFAULT_TO
) {
  query <- paste0(
    "orgId=1",
    "&panelId=", GRAFANA_SALINITY_PANEL_ID,
    "&var-metric=mean",
    "&var-Buoy=", URLencode(buoy_label, reserved = TRUE),
    "&from=", from,
    "&to=", to
  )
  paste0(
    GRAFANA_BASE,
    "/d-solo/",
    GRAFANA_DASHBOARD_UID,
    "/",
    GRAFANA_DASHBOARD_SLUG,
    "?",
    query
  )
}

grafana_dashboard_url <- function(
    from = GRAFANA_DEFAULT_FROM,
    to = GRAFANA_DEFAULT_TO
) {
  paste0(
    GRAFANA_BASE,
    "/d/",
    GRAFANA_DASHBOARD_UID,
    "/",
    GRAFANA_DASHBOARD_SLUG,
    "?orgId=1&from=",
    from,
    "&to=",
    to,
    "&var-metric=mean&viewPanel=panel-",
    GRAFANA_SALINITY_PANEL_ID
  )
}

grafana_salinity_query_body <- function(
    buoy_label = NULL,
    from = GRAFANA_DEFAULT_FROM,
    to = GRAFANA_DEFAULT_TO,
    interval_ms = 86400000L,
    max_data_points = 1000L
) {
  if (is.null(buoy_label)) {
    query <- paste(
      "SELECT mean(\"sea_water_practical_salinity\")",
      "FROM \"ndbc_buoy\"",
      "WHERE $timeFilter",
      "GROUP BY time($__interval), \"location\" fill(null)"
    )
  } else {
    query <- sprintf(
      paste(
        "SELECT mean(\"sea_water_practical_salinity\")",
        "FROM \"ndbc_buoy\"",
        "WHERE (\"location\" =~ /^%s$/) AND $timeFilter",
        "GROUP BY time($__interval), \"location\" fill(null)"
      ),
      grafana_regex_buoy(buoy_label)
    )
  }

  list(
    queries = list(
      list(
        refId = "A",
        datasource = list(
          type = "influxdb",
          uid = GRAFANA_INFLUX_DATASOURCE_UID
        ),
        query = query,
        rawQuery = TRUE,
        intervalMs = interval_ms,
        maxDataPoints = max_data_points
      )
    ),
    from = as.character(from),
    to = as.character(to)
  )
}

grafana_parse_timeseries_frame <- function(frame) {
  if (is.null(frame$data) || is.null(frame$data$values)) {
    return(data.frame(
      time = as.POSIXct(character()),
      salinity = numeric(),
      location = character(),
      stringsAsFactors = FALSE
    ))
  }

  values <- frame$data$values
  if (length(values) < 2 || length(values[[1]]) == 0) {
    return(data.frame(
      time = as.POSIXct(character()),
      salinity = numeric(),
      location = character(),
      stringsAsFactors = FALSE
    ))
  }

  location <- frame$schema$fields[[2]]$labels$location
  if (is.null(location)) {
    location <- NA_character_
  }

  time_ms <- suppressWarnings(as.numeric(unlist(values[[1]], use.names = FALSE)))
  salinity <- suppressWarnings(as.numeric(unlist(values[[2]], use.names = FALSE)))
  valid <- !is.na(time_ms) & !is.na(salinity)
  if (!any(valid)) {
    return(data.frame(
      time = as.POSIXct(character()),
      salinity = numeric(),
      location = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    time = as.POSIXct(time_ms[valid] / 1000, origin = "1970-01-01", tz = "UTC"),
    salinity = as.numeric(salinity[valid]),
    location = location,
    stringsAsFactors = FALSE
  )
}

fetch_grafana_salinity_timeseries <- function(
    buoy_label = NULL,
    from = GRAFANA_DEFAULT_FROM,
    to = GRAFANA_DEFAULT_TO
) {
  body <- grafana_salinity_query_body(
    buoy_label = buoy_label,
    from = from,
    to = to
  )

  resp <- httr2::request(paste0(GRAFANA_BASE, "/api/ds/query")) |>
    httr2::req_method("POST") |>
    httr2::req_body_json(body) |>
    httr2::req_perform()

  payload <- httr2::resp_body_json(resp)
  frames <- payload$results$A$frames
  if (is.null(frames) || length(frames) == 0) {
    return(data.frame(
      time = as.POSIXct(character()),
      salinity = numeric(),
      location = character(),
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(buoy_label)) {
    dplyr::bind_rows(lapply(frames, grafana_parse_timeseries_frame))
  } else {
    grafana_parse_timeseries_frame(frames[[1]])
  }
}

plot_grafana_salinity_panel <- function(
    data,
    buoy_label,
    solo_url,
    low_threshold = 25,
    high_threshold = 35
) {
  if (nrow(data) == 0) {
    return(htmltools::tags$div(
      class = "grafana-ts-card",
      htmltools::tags$div(class = "grafana-ts-title", buoy_label),
      htmltools::tags$p(class = "grafana-ts-empty", "No salinity data returned for this buoy."),
      htmltools::tags$a(
        href = solo_url,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open Grafana panel"
      )
    ))
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = time, y = salinity)
  ) +
    ggplot2::geom_hline(
      yintercept = low_threshold,
      linetype = "dashed",
      color = "#2166ac",
      linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = high_threshold,
      linetype = "dashed",
      color = "#b2182b",
      linewidth = 0.4
    ) +
    ggplot2::geom_line(color = "#1f6feb", linewidth = 0.7) +
    ggplot2::scale_x_datetime(date_labels = "%Y-%m") +
    ggplot2::labs(
      x = NULL,
      y = "Salinity (PSU)",
      title = sprintf("In-situ salinity at ~1 m depth from %s", buoy_label)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  htmltools::tags$div(
    class = "grafana-ts-card",
    htmltools::tags$div(
      class = "grafana-ts-links",
      htmltools::tags$a(
        href = solo_url,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open Grafana panel"
      )
    ),
    htmltools::plotTag(p, alt = sprintf("Salinity time series for %s", buoy_label))
  )
}

plot_grafana_salinity_panels <- function(stations, salinity_ts, from, to) {
  cards <- lapply(stations, function(station_info) {
    buoy_label <- grafana_buoy_label(station_info)
    station_data <- salinity_ts |>
      dplyr::filter(.data$location == buoy_label)
    solo_url <- grafana_solo_panel_url(buoy_label, from = from, to = to)
    plot_grafana_salinity_panel(station_data, buoy_label, solo_url)
  })

  htmltools::tags$div(class = "grafana-ts-grid", cards)
}
