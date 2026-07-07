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

field_raster <- function(field) {
  raster::rasterFromXYZ(
    field[, c("longitude", "latitude", "salinity")],
    crs = "+proj=longlat +datum=WGS84"
  )
}

plot_salinity_field <- function(field, observations) {
  if (is.null(field) || nrow(field) == 0) {
    return(NULL)
  }

  obs <- as.data.frame(observations)
  salinity_range <- range(
    c(field$salinity, obs$salinity),
    na.rm = TRUE
  )
  pal <- leaflet::colorNumeric(
    viridisLite::turbo(256),
    domain = salinity_range,
    na.color = "transparent"
  )

  lon_min <- min(field$longitude, na.rm = TRUE)
  lon_max <- max(field$longitude, na.rm = TRUE)
  lat_min <- min(field$latitude, na.rm = TRUE)
  lat_max <- max(field$latitude, na.rm = TRUE)

  field_rast <- field_raster(field)

  map <- leaflet::leaflet(
    width = "100%",
    height = 520
  ) |>
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
    leaflet::addCircleMarkers(
      data = obs,
      lng = ~longitude,
      lat = ~latitude,
      radius = 6,
      fillColor = ~pal(salinity),
      fillOpacity = 0.95,
      color = "#1a1a1a",
      weight = 1,
      stroke = TRUE,
      label = ~sprintf(
        "%s: %.1f PSU\n%s",
        location,
        salinity,
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
      values = salinity_range,
      title = "Salinity (PSU)",
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
