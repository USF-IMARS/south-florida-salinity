SALINITY_GAUGE_MIN <- 10
SALINITY_GAUGE_MAX <- 50
SALINITY_LOW_THRESHOLD <- 25
SALINITY_HIGH_THRESHOLD <- 35

salinity_gauge_status <- function(value) {
  if (is.na(value)) {
    return("missing")
  }
  if (value < SALINITY_LOW_THRESHOLD) {
    return("low")
  }
  if (value > SALINITY_HIGH_THRESHOLD) {
    return("high")
  }
  "normal"
}

salinity_gauge_status_label <- function(status) {
  switch(
    status,
    low = "Low",
    normal = "Normal",
    high = "High",
    missing = "No data"
  )
}

salinity_gauge_position_pct <- function(value) {
  if (is.na(value)) {
    return(NA_real_)
  }
  pct <- (value - SALINITY_GAUGE_MIN) / (SALINITY_GAUGE_MAX - SALINITY_GAUGE_MIN) * 100
  max(0, min(100, pct))
}

salinity_threshold_pct <- function(threshold) {
  (threshold - SALINITY_GAUGE_MIN) / (SALINITY_GAUGE_MAX - SALINITY_GAUGE_MIN) * 100
}

salinity_gauge <- function(
    station,
    location,
    value,
    time,
    station_id = NA_character_
) {
  status <- salinity_gauge_status(value)
  position <- salinity_gauge_position_pct(value)
  low_pct <- salinity_threshold_pct(SALINITY_LOW_THRESHOLD)
  high_pct <- salinity_threshold_pct(SALINITY_HIGH_THRESHOLD)

  value_label <- if (is.na(value)) {
    "--"
  } else {
    sprintf("%.1f PSU", value)
  }

  time_label <- if (inherits(time, "POSIXct") || inherits(time, "POSIXt")) {
    format(time, "%Y-%m-%d %H:%M UTC")
  } else if (is.na(time)) {
    ""
  } else {
    as.character(time)
  }

  marker <- if (is.na(position)) {
    htmltools::tags$div(class = "salinity-gauge-marker salinity-gauge-marker-missing")
  } else {
    htmltools::tags$div(
      class = "salinity-gauge-marker",
      style = sprintf("left: %.2f%%;", position),
      title = sprintf("%.1f PSU", value)
    )
  }

  htmltools::tags$div(
    class = "salinity-gauge-card",
    htmltools::tags$div(class = "salinity-gauge-title", station),
    htmltools::tags$div(
      class = "salinity-gauge-subtitle",
      location
    ),
    htmltools::tags$div(
      class = "salinity-gauge-track",
      style = sprintf(
        "--low-end:%.2f%%; --high-start:%.2f%%;",
        low_pct,
        high_pct
      ),
      htmltools::tags$div(class = "salinity-gauge-zone salinity-gauge-zone-low"),
      htmltools::tags$div(class = "salinity-gauge-zone salinity-gauge-zone-normal"),
      htmltools::tags$div(class = "salinity-gauge-zone salinity-gauge-zone-high"),
      htmltools::tags$div(class = "salinity-gauge-threshold salinity-gauge-threshold-low"),
      htmltools::tags$div(class = "salinity-gauge-threshold salinity-gauge-threshold-high"),
      marker
    ),
    htmltools::tags$div(
      class = "salinity-gauge-scale",
      htmltools::tags$span(sprintf("%d", SALINITY_GAUGE_MIN)),
      htmltools::tags$span(sprintf("%d", SALINITY_LOW_THRESHOLD)),
      htmltools::tags$span(sprintf("%d", SALINITY_HIGH_THRESHOLD)),
      htmltools::tags$span(sprintf("%d", SALINITY_GAUGE_MAX))
    ),
    htmltools::tags$div(class = "salinity-gauge-value", value_label),
    htmltools::tags$div(
      class = sprintf("salinity-gauge-status salinity-gauge-status-%s", status),
      salinity_gauge_status_label(status)
    ),
    htmltools::tags$div(class = "salinity-gauge-time", time_label)
  )
}

plot_salinity_gauges <- function(data) {
  cards <- lapply(seq_len(nrow(data)), function(i) {
    row <- data[i, ]
    salinity_gauge(
      station = row$station,
      location = row$location,
      value = row$salinity,
      time = row$time,
      station_id = row$station_id
    )
  })

  htmltools::tags$div(class = "salinity-gauge-grid", cards)
}

salinity_gauge_legend <- function() {
  htmltools::tags$div(
    class = "salinity-gauge-legend",
    htmltools::tags$div(
      class = "salinity-gauge-legend-item",
      htmltools::tags$span(class = "salinity-gauge-legend-swatch salinity-gauge-zone-low"),
      sprintf("Low: < %.0f PSU", SALINITY_LOW_THRESHOLD)
    ),
    htmltools::tags$div(
      class = "salinity-gauge-legend-item",
      htmltools::tags$span(class = "salinity-gauge-legend-swatch salinity-gauge-zone-normal"),
      sprintf(
        "Normal: %.0f–%.0f PSU",
        SALINITY_LOW_THRESHOLD,
        SALINITY_HIGH_THRESHOLD
      )
    ),
    htmltools::tags$div(
      class = "salinity-gauge-legend-item",
      htmltools::tags$span(class = "salinity-gauge-legend-swatch salinity-gauge-zone-high"),
      sprintf("High: > %.0f PSU", SALINITY_HIGH_THRESHOLD)
    )
  )
}
