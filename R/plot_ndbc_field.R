field_slice_grid <- function(field_slice) {
  lons <- sort(unique(field_slice$longitude))
  lats <- sort(unique(field_slice$latitude))
  n_lon <- length(lons)
  n_lat <- length(lats)

  dx <- if (n_lon > 1) diff(lons)[1] / 2 else 0.01
  dy <- if (n_lat > 1) diff(lats)[1] / 2 else 0.01

  field_slice$lon1 <- field_slice$longitude - dx
  field_slice$lon2 <- field_slice$longitude + dx
  field_slice$lat1 <- field_slice$latitude - dy
  field_slice$lat2 <- field_slice$latitude + dy
  field_slice
}

field_raster <- function(field, value_col) {
  xyz <- field[, c("longitude", "latitude", value_col)]
  names(xyz)[3] <- "value"
  raster::rasterFromXYZ(
    xyz,
    crs = "+proj=longlat +datum=WGS84"
  )
}

plot_ndbc_field <- function(
    field,
    observations,
    value_col,
    legend_title,
    popup_suffix = "",
    palette = viridisLite::turbo(256)
) {
  if (is.null(field) || nrow(field) == 0) {
    return(NULL)
  }

  obs <- as.data.frame(observations)
  obs$marker_value <- obs[[value_col]]
  value_range <- range(
    c(field[[value_col]], obs$marker_value),
    na.rm = TRUE
  )
  pal <- leaflet::colorNumeric(
    palette,
    domain = value_range,
    na.color = "transparent"
  )

  lon_min <- min(field$longitude, na.rm = TRUE)
  lon_max <- max(field$longitude, na.rm = TRUE)
  lat_min <- min(field$latitude, na.rm = TRUE)
  lat_max <- max(field$latitude, na.rm = TRUE)

  field_rast <- field_raster(field, value_col)
  grid <- field_slice_grid(field)
  grid$popup_label <- sprintf(
    paste0("%.1f", popup_suffix),
    grid[[value_col]]
  )

  map <- leaflet::leaflet(
    width = "100%",
    height = 520,
    options = leaflet::leafletOptions(preferCanvas = TRUE)
  ) |>
    leaflet::addMapPane("fieldValues", zIndex = 620) |>
    leaflet::addProviderTiles(
      leaflet::providers$Esri.OceanBasemap,
      group = "Ocean"
    ) |>
    leaflet::addProviderTiles(
      leaflet::providers$CartoDB.Positron,
      group = "Light"
    ) |>
    leaflet::addProviderTiles(
      leaflet::providers$OpenStreetMap,
      group = "OpenStreetMap"
    ) |>
    leaflet::addRasterImage(
      field_rast,
      colors = pal,
      opacity = 0.85,
      group = "Interpolated field"
    ) |>
    leaflet::addRectangles(
      data = grid,
      lng1 = ~lon1,
      lat1 = ~lat1,
      lng2 = ~lon2,
      lat2 = ~lat2,
      fillColor = "#000000",
      fillOpacity = 0.001,
      color = "#000000",
      opacity = 0,
      weight = 0,
      popup = ~popup_label,
      group = "Interpolated field",
      options = leaflet::pathOptions(
        pane = "fieldValues",
        interactive = TRUE
      ),
      highlightOptions = leaflet::highlightOptions(
        weight = 1,
        color = "#333333",
        fillOpacity = 0.15,
        bringToFront = TRUE
      )
    ) |>
    leaflet::addCircleMarkers(
      data = obs,
      lng = ~longitude,
      lat = ~latitude,
      radius = 6,
      fillColor = ~pal(marker_value),
      fillOpacity = 0.95,
      color = "#1a1a1a",
      weight = 1,
      stroke = TRUE,
      label = ~sprintf(
        "%s: %.1f%s\n%s",
        location,
        marker_value,
        popup_suffix,
        format(time, "%Y-%m-%d %H:%M UTC")
      ),
      group = "Observations"
    ) |>
    leaflet::addLayersControl(
      baseGroups = c("Ocean", "Light", "OpenStreetMap"),
      overlayGroups = c("Interpolated field", "Observations"),
      options = leaflet::layersControlOptions(collapsed = TRUE)
    ) |>
    leaflet::addLegend(
      pal = pal,
      values = value_range,
      title = legend_title,
      position = "bottomright"
    ) |>
    leaflet::addScaleBar(position = "bottomleft") |>
    leaflet::fitBounds(
      lng1 = lon_min,
      lat1 = lat_min,
      lng2 = lon_max,
      lat2 = lat_max
    )

  map
}

plot_salinity_field <- function(field, observations) {
  plot_ndbc_field(
    field,
    observations,
    value_col = "salinity",
    legend_title = "Salinity (PSU)",
    popup_suffix = " PSU"
  )
}
