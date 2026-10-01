# ============================================================
# GUIDED FROM ABOVE  (v2)
# GeoHipster 2027 calendar submission
#
# GPS adoption on U.S. cropland (USDA ERS EIB-248), 2013 -> 2019,
# drawn as automated-guidance passes under the GPS constellation
# over Manhattan, Kansas during the May 2024 Gannon storm.
#
# Encoding
#   * Each state is "planted" with parallel N-S guidance passes.
#     Pass density rises with the 2019 share of cropland acres
#     using GPS, rescaled over the observed range (denser = higher;
#     NOT proportional -- the legend keys show real data values).
#   * Pass colour = percentage-point change in that share,
#     2013-2019, capped at +/- `cap` (90th pct of |change|).
#   * State outline = the headland pass.
#   * Sky = GPS satellites above Manhattan, KS at `obs_time`
#     (panorama: azimuth left->right, elevation bottom->top).
#
# Output: 10.5 x 8 in landscape, 600 dpi PNG (GeoHipster spec).
# Everything important sits >= 0.75 in from the edge
# (0.25 in bleed + 0.5 in margin). Bottom-left corner is left
# clear for the GeoHipster author/logo box.
# ============================================================

library(tidyverse)   # dplyr, ggplot2, readr, stringr, purrr, lubridate
library(sf)
library(tigris)

options(tigris_use_cache = TRUE)
set.seed(2027)

# ============================================================
# 0. CANVAS, COLOURS, FONTS
# ============================================================

canvas_w_in <- 10.5
canvas_h_in <- 8.0
export_dpi  <- 600
safe_in     <- 0.75            # bleed + margin
out_file    <- "guided_from_above.png"

crs_conus <- 5070              # NAD83 / CONUS Albers Equal Area

col_bg     <- "#05070F"        # night
col_land   <- "#0A1122"        # unplowed ground
col_edge   <- "#22304F"        # outlines of states with no estimate
col_grat   <- "#0F1830"        # graticule
col_star   <- "#DDE6FF"
col_orbit  <- "#3A4E7A"
col_signal <- "#7FD8FF"
col_text   <- "#E8ECF7"
col_sub    <- "#8A96B8"
col_na     <- "#4A5570"

# Diverging ramp on dark ground: loss -> no change -> gain
ramp_cols <- c("#FF4F9A", "#7FD8FF", "#FFE08A", "#FFB000")
ramp_vals <- c(0, 0.5, 0.78, 1)

# Fonts: Space Grotesk / Space Mono via showtext if available,
# otherwise fall back to system sans / mono.
font_title <- "sans"
font_mono  <- "mono"
use_showtext <- FALSE
if (requireNamespace("showtext", quietly = TRUE)) {
  use_showtext <- tryCatch({
    sysfonts::font_add_google("Space Grotesk", "spacegrotesk")
    sysfonts::font_add_google("Space Mono", "spacemono")
    TRUE
  }, error = function(e) FALSE)
}
if (use_showtext) {
  showtext::showtext_auto()
  # Text sizes are calibrated for the saved file; the RStudio
  # preview will look oversized. Judge the PNG, not the viewer.
  showtext::showtext_opts(dpi = export_dpi)
  font_title <- "spacegrotesk"
  font_mono  <- "spacemono"
}

pt <- function(x) x / .pt      # points -> ggplot text size (mm)

# ============================================================
# 1. DATA (USDA ERS state estimates)
# ============================================================

read_ers <- function(path, yr) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      year = yr,
      fips = str_pad(as.character(fips), width = 2, side = "left", pad = "0")
    )
}

df_2013 <- read_ers("data/fig3_data_2013.csv", 2013)
df_2019 <- read_ers("data/fig3_data_2019.csv", 2019)

# full_join so a state with a 2019 estimate but no 2013 estimate
# still gets passes (drawn in neutral grey, change = NA)
df_diff <- df_2019 %>%
  select(fips, region, shareacres_2019 = shareacres, sharefarms_2019 = sharefarms) %>%
  full_join(
    df_2013 %>% select(fips, shareacres_2013 = shareacres),
    by = "fips"
  ) %>%
  mutate(change_acres = shareacres_2019 - shareacres_2013)

# ============================================================
# 2. GEOMETRY
# ============================================================

states <- tigris::states(cb = TRUE, year = 2024) %>%
  filter(STUSPS %in% state.abb, !STUSPS %in% c("AK", "HI")) %>%
  st_transform(crs_conus) %>%
  select(GEOID, STUSPS, NAME)

map_sf <- states %>%
  left_join(df_diff, by = c("GEOID" = "fips"))

# Canvas in projected metres. The map is shifted left so the
# Atlantic below New England holds the legend and title block.
bb    <- st_bbox(states)
map_w <- unname(bb["xmax"] - bb["xmin"])
upi   <- map_w * 1.40 / canvas_w_in        # metres per printed inch
cx0   <- unname(bb["xmin"]) - 0.90 * upi
cx1   <- cx0 + canvas_w_in * upi
cy0   <- unname(bb["ymin"]) - 0.80 * upi
cy1   <- cy0 + canvas_h_in * upi
safe  <- safe_in * upi
inch  <- function(x) x * upi               # inches -> metres

horizon_y <- unname(bb["ymax"]) + inch(0.12)
sky_top   <- cy1 - safe

# ============================================================
# 3. GUIDANCE PASSES
# ============================================================

has_level   <- map_sf %>% filter(!is.na(shareacres_2019))
no_estimate <- map_sf %>% filter(is.na(shareacres_2019))

# Pass spacing: lines per printed inch rescaled linearly over the
# observed 2019 range -- lowest state 1 pass / 0.30 in, highest
# 1 pass / 0.03 in. Defined HERE (before the passes are built) so
# the map and the legend key use the same function.
sh_rng <- range(has_level$shareacres_2019, na.rm = TRUE)
pass_spacing <- function(share) {
  lpi <- scales::rescale(share, to = c(1 / 0.30, 1 / 0.03), from = sh_rng)
  inch(1 / lpi)
}

make_passes <- function(geom, spacing) {
  b    <- st_bbox(geom)
  xmin <- b[["xmin"]]; xmax <- b[["xmax"]]
  # snap to a global grid so neighbouring states line up
  x0 <- ceiling(xmin / spacing) * spacing
  xs <- if (x0 <= xmax) seq(x0, xmax, by = spacing) else numeric(0)
  # State narrower than one pass spacing: give it a single pass
  # through its interior so a low-but-nonzero share isn't drawn
  # the same as "no estimate"
  if (length(xs) == 0) {
    xs <- st_coordinates(st_point_on_surface(geom))[1, "X"]
  }
  ml <- st_multilinestring(lapply(xs, function(x)
    matrix(c(x, x, b[["ymin"]] - 1, b[["ymax"]] + 1), ncol = 2)))
  out <- st_intersection(st_sfc(ml, crs = st_crs(geom)), geom)
  if (length(out) == 0) return(st_sfc(st_multilinestring(), crs = st_crs(geom)))
  if (inherits(out, "sfc_GEOMETRYCOLLECTION")) out <- st_collection_extract(out, "LINESTRING")
  st_union(st_cast(out, "MULTILINESTRING"))
}

# Colour: cap at the 90th percentile of |change| (rounded up to 5 pp)
cap <- unname(ceiling(quantile(abs(map_sf$change_acres), 0.9, na.rm = TRUE) / 5) * 5)
change_rng <- c(-cap, cap)
change_col <- function(x) {
  r   <- scales::rescale(pmin(pmax(x, -cap), cap), to = c(0, 1), from = change_rng)
  out <- scales::gradient_n_pal(ramp_cols, ramp_vals)(r)
  out[is.na(x)] <- col_na
  out
}

# States beyond the cap (worth checking against ERS reliability flags)
map_sf %>% st_drop_geometry() %>%
  filter(abs(change_acres) > cap) %>%
  select(STUSPS, shareacres_2013, shareacres_2019, change_acres) %>%
  print()

passes <- has_level %>%
  mutate(
    colour   = change_col(change_acres),
    geometry = do.call(c, map2(st_geometry(.), shareacres_2019,
                               ~ make_passes(st_sfc(.x, crs = crs_conus), pass_spacing(.y))))
  ) %>%
  st_set_geometry("geometry") %>%
  select(GEOID, colour, geometry)

headlands <- has_level %>%
  mutate(colour = change_col(change_acres)) %>%
  st_boundary()

graticule <- st_graticule(
  st_bbox(st_transform(states, 4326)),
  crs = st_crs(crs_conus), lon = seq(-125, -65, 5), lat = seq(25, 50, 5)
)

# ============================================================
# 4. GPS CONSTELLATION ABOVE MANHATTAN, KS
# ============================================================

# Gannon storm: G5 conditions from ~17 UTC 10 May 2024. 22:00 UTC is
# 5 pm CDT (daylight, no stars). 04:00 UTC 11 May = 11 pm CDT 10 May:
# dark, storm ongoing, aurora reported across Kansas.
obs_time <- as.POSIXct("2024-05-11 04:00:00", tz = "UTC")
obs_lat  <- 39.1836
obs_lon  <- -96.5717
obs_h_km <- 0.32
elev_mask <- 15                                             # deg, for signal lines

# Historical element sets for the storm week, downloaded from
# Space-Track.org (class gp_history, OBJECT_NAME ~ "GPS",
# EPOCH 2024-05-07--2024-05-14, format CSV). Same OMM columns as
# CelesTrak. If this file is missing the script falls back to
# current CelesTrak elements and warns that the sky is approximate.
history_file <- "data/gps_gp_history_2024-05.csv"

# Optional: PRNs set healthy on the obs date (from NAVCEN / IGS
# status for May 2024). NULL = keep every GPS object in the file.
active_prns <- NULL

mu_E <- 398600.4418; R_E <- 6378.137; J2 <- 1.08262668e-3

nominal_constellation <- function(t) {
  # 24-slot baseline (6 planes x 4); used only if no elements are available
  expand_grid(plane = 1:6, slot = 1:4) %>%
    transmute(
      label = paste0(LETTERS[plane], slot),
      epoch = t, n_revday = 2.00563, e = 0, i = 55,
      raan = (plane - 1) * 60, argp = 0,
      M = ((slot - 1) * 90 + (plane - 1) * 15) %% 360,
      source = "nominal"
    )
}

tidy_omm <- function(el, source) {
  el %>%
    transmute(
      prn   = str_match(OBJECT_NAME, "PRN\\s*(\\d+)")[, 2],
      label = if_else(is.na(prn), OBJECT_NAME, paste0("G", str_pad(prn, 2, pad = "0"))),
      epoch = lubridate::ymd_hms(EPOCH, tz = "UTC"),
      n_revday = MEAN_MOTION, e = ECCENTRICITY, i = INCLINATION,
      raan = RA_OF_ASC_NODE, argp = ARG_OF_PERICENTER, M = MEAN_ANOMALY,
      source = source
    )
}

get_gps_elements <- function(t) {
  omm_cols <- cols(EPOCH = col_character(), .default = col_guess())
  if (file.exists(history_file)) {
    el <- read_csv(history_file, col_types = omm_cols) %>%
      mutate(epoch_t = lubridate::ymd_hms(EPOCH, tz = "UTC")) %>%
      group_by(NORAD_CAT_ID) %>%
      slice_min(abs(as.numeric(difftime(epoch_t, t, units = "secs"))),
                n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      tidy_omm("Space-Track")
    # NB: OBJECT_NAME carries the catalogue's *current* PRN; check any
    # vehicle whose PRN was reassigned after May 2024.
    if (!is.null(active_prns)) {
      el <- el %>% filter(label %in% paste0("G", str_pad(active_prns, 2, pad = "0")))
    }
  } else {
    el <- tryCatch(
      read_csv("https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=csv",
               col_types = omm_cols) %>% tidy_omm("CelesTrak"),
      error = function(e) NULL
    )
    if (is.null(el) || nrow(el) < 20) {
      message("No element sets available - using nominal constellation")
      return(nominal_constellation(t))
    }
  }
  gap_days <- max(abs(as.numeric(difftime(el$epoch, t, units = "days"))))
  if (gap_days > 30) {
    warning(sprintf(
      "Element epochs are up to %.0f days from obs_time: the sky is NOT the real sky at that moment.",
      gap_days), call. = FALSE)
  }
  el %>% mutate(gap_days = gap_days)
}

solve_kepler <- function(M, e) {
  E <- M
  for (k in 1:12) E <- E - (E - e * sin(E) - M) / (1 - e * cos(E))
  E
}

gmst_rad <- function(t) {
  jd <- as.numeric(t) / 86400 + 2440587.5
  ((280.46061837 + 360.98564736629 * (jd - 2451545.0)) %% 360) * pi / 180
}

observer_ecef <- function(lat, lon, h) {
  f <- 1 / 298.257223563; e2 <- f * (2 - f)
  la <- lat * pi / 180; lo <- lon * pi / 180
  N <- R_E / sqrt(1 - e2 * sin(la)^2)
  c((N + h) * cos(la) * cos(lo), (N + h) * cos(la) * sin(lo), (N * (1 - e2) + h) * sin(la))
}

# Two-body propagation with secular J2 drift of RAAN and argument
# of perigee. Good to roughly a degree within a few weeks of the
# element epoch - plenty for a sky picture, not for navigation.
sat_azel <- function(el, times) {
  o  <- observer_ecef(obs_lat, obs_lon, obs_h_km)
  la <- obs_lat * pi / 180; lo <- obs_lon * pi / 180
  expand_grid(el, time = times) %>%
    mutate(
      dt   = as.numeric(difftime(time, epoch, units = "secs")),
      n    = n_revday * 2 * pi / 86400,
      a    = (mu_E / n^2)^(1 / 3),
      p    = a * (1 - e^2),
      ir   = i * pi / 180,
      Om   = raan * pi / 180 - 1.5 * n * J2 * (R_E / p)^2 * cos(ir) * dt,
      w    = argp * pi / 180 + 0.75 * n * J2 * (R_E / p)^2 * (5 * cos(ir)^2 - 1) * dt,
      Mt   = (M * pi / 180 + n * dt) %% (2 * pi),
      E    = solve_kepler(Mt, e),
      nu   = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2)),
      r    = a * (1 - e * cos(E)),
      u    = w + nu,
      xi   = r * (cos(Om) * cos(u) - sin(Om) * sin(u) * cos(ir)),
      yi   = r * (sin(Om) * cos(u) + cos(Om) * sin(u) * cos(ir)),
      zi   = r * sin(u) * sin(ir),
      g    = gmst_rad(time),
      dx   =  cos(g) * xi + sin(g) * yi - o[1],
      dy   = -sin(g) * xi + cos(g) * yi - o[2],
      dz   =  zi - o[3],
      east  = -sin(lo) * dx + cos(lo) * dy,
      north = -sin(la) * cos(lo) * dx - sin(la) * sin(lo) * dy + cos(la) * dz,
      up    =  cos(la) * cos(lo) * dx + cos(la) * sin(lo) * dy + sin(la) * dz,
      az   = (atan2(east, north) * 180 / pi) %% 360,
      elev = asin(up / sqrt(dx^2 + dy^2 + dz^2)) * 180 / pi
    ) %>%
    select(label, source, time, az, elev)
}

gps_el     <- get_gps_elements(obs_time)
sky_source <- unique(gps_el$source)
# The sky is "real" only with element sets close to obs_time
sky_is_real <- sky_source != "nominal" && isTRUE(all(gps_el$gap_days <= 30))

# Panorama: azimuth 0-360 across the safe width (S in the middle),
# elevation 0-90 from the horizon to the top of the safe area.
sky_x <- function(az)   (cx0 + safe) + az / 360 * (cx1 - cx0 - 2 * safe)
sky_y <- function(elev) horizon_y + inch(0.10) + elev / 90 * (sky_top - horizon_y - inch(0.10))

sats_now <- sat_azel(gps_el, obs_time) %>%
  filter(elev > 0) %>%
  mutate(x = sky_x(az), y = sky_y(elev))

tracks <- sat_azel(gps_el, seq(obs_time - 3 * 3600, obs_time + 3 * 3600, by = 300)) %>%
  arrange(label, time) %>%
  group_by(label) %>%
  mutate(
    brk = elev <= 0 | abs(az - lag(az, default = first(az))) > 180,
    seg = cumsum(brk)
  ) %>%
  filter(elev > 0) %>%
  ungroup() %>%
  mutate(x = sky_x(az), y = sky_y(elev), grp = paste(label, seg))

home <- st_sfc(st_point(c(obs_lon, obs_lat)), crs = 4326) %>%
  st_transform(crs_conus) %>% st_coordinates()

signals <- sats_now %>%
  filter(elev >= elev_mask) %>%
  mutate(xend = home[1], yend = home[2])

message(sprintf("%d satellites above horizon, %d above %d deg mask (%s)",
                nrow(sats_now), nrow(signals), elev_mask, sky_source))

# ============================================================
# 5. SKY DRESSING
# ============================================================

stars <- tibble(
  x = runif(1100, cx0, cx1),
  y = runif(1100, horizon_y, cy1),
  size  = rexp(1100, 9) + 0.02,
  alpha = runif(1100, 0.15, 0.9)
)

horizon_glow <- tibble(k = 0:39) %>%
  mutate(
    ymin  = horizon_y + inch(k * 0.02),
    ymax  = ymin + inch(0.02),
    alpha = 0.22 * exp(-k / 9)
  )

# Aurora: decorative glow on the northern horizon (both ends of the
# panorama), green low, magenta high. Set FALSE to drop it.
draw_aurora <- TRUE
aurora <- expand_grid(az = seq(-60, 60, by = 1), k = 0:44) %>%
  mutate(
    strength = exp(-(az / 30)^2) * (0.75 + 0.25 * sin(az / 5)),  # curtain ripple
    y0   = horizon_y + inch(k * 0.025),
    y1   = y0 + inch(0.025),
    fill = scales::gradient_n_pal(c("#35F2A0", "#35F2A0", "#B45CFF", "#FF5FD2"),
                                  c(0, 0.35, 0.75, 1))(k / 44),
    alpha = 0.16 * strength * exp(-k / 18),
    az360 = az %% 360
  ) %>%
  # az -60..60 wraps to 300..360 (right edge) and 0..60 (left edge)
  mutate(xmin = sky_x(az360), xmax = xmin + (cx1 - cx0 - 2 * safe) / 360 * 1.3)

compass <- tibble(
  az  = c(0, 90, 180, 270, 360),
  lab = c("N", "E", "S", "W", "N")
) %>% mutate(x = sky_x(az))

# ============================================================
# 6. LEGEND BLOCK (bottom right, inside the safe area)
# ============================================================

lg_w  <- inch(1.90)
lg_x1 <- cx1 - safe
lg_x0 <- lg_x1 - lg_w
lg_y0 <- cy0 + safe

# 6a. colour ramp for change
ramp_y  <- lg_y0 + inch(0.52)
ramp_h  <- inch(0.09)
ramp_df <- tibble(v = seq(change_rng[1], change_rng[2], length.out = 120)) %>%
  mutate(
    xmin = lg_x0 + (row_number() - 1) / n() * lg_w,
    xmax = xmin + 1.6 * lg_w / n(),   # overlap to hide seams
    fill = change_col(v)
  )
# Explicit breaks so the capped ends (<= / >=) are always labelled
ramp_breaks <- c(-cap, -cap / 2, 0, cap / 2, cap)
fmt <- function(v) format(v, trim = TRUE, drop0trailing = TRUE)
ramp_ticks <- tibble(
  v = ramp_breaks,
  x = lg_x0 + (v - change_rng[1]) / diff(change_rng) * lg_w,
  lab = case_when(v <= -cap ~ paste0("\u2264", fmt(v)),
                  v >=  cap ~ paste0("\u2265+", fmt(v)),
                  v > 0     ~ paste0("+", fmt(v)),
                  TRUE      ~ fmt(v))
)

# 6b. pass-spacing key: same pass_spacing() as the map, keys at the
#     10th / 50th / 90th percentile of the actual 2019 shares
key_side  <- inch(0.42)
key_gap   <- (lg_w - 3 * key_side) / 2
key_y0    <- lg_y0 + inch(1.10)
key_share <- unname(round(quantile(has_level$shareacres_2019, c(0.1, 0.5, 0.9))))

key_boxes <- map(0:2, function(k) {
  x0 <- lg_x0 + k * (key_side + key_gap)
  st_polygon(list(matrix(c(x0, key_y0, x0 + key_side, key_y0,
                           x0 + key_side, key_y0 + key_side,
                           x0, key_y0 + key_side, x0, key_y0),
                         ncol = 2, byrow = TRUE)))
}) %>% st_sfc(crs = crs_conus)
key_passes <- do.call(c, map2(seq_along(key_boxes), key_share,
                              ~ make_passes(key_boxes[.x], pass_spacing(.y))))
key_labels <- tibble(
  x = map_dbl(seq_along(key_boxes), ~ st_bbox(key_boxes[.x])[["xmin"]]) + key_side / 2,
  lab = paste0(key_share, "%")
)

sky_credit <- switch(sky_source, "Space-Track" = "Space-Track.org",
                     "CelesTrak" = "CelesTrak", "nominal" = "nominal GPS constellation")
credit <- paste0("USDA ERS EIB-248 (McFadden et al. 2023)\n",
                 "Census TIGER \u00b7 ", sky_credit, " \u00b7 @SpacePlowboy")

# ============================================================
# 7. COMPOSE
# ============================================================

pass_layers <- function(data) {
  list(
    geom_sf(data = data, aes(colour = colour), linewidth = 0.95, alpha = 0.07),
    geom_sf(data = data, aes(colour = colour), linewidth = 0.45, alpha = 0.22),
    geom_sf(data = data, aes(colour = colour), linewidth = 0.16, alpha = 1)
  )
}

txt <- function(x, y, label, size, family = font_mono, colour = col_sub,
                hjust = 0, vjust = 0, fontface = "plain", lineheight = 0.95) {
  annotate("text", x = x, y = y, label = label, size = pt(size), family = family,
           colour = colour, hjust = hjust, vjust = vjust, fontface = fontface,
           lineheight = lineheight)
}

# Satellite labels: ggrepel if installed (avoids G12/G15-type collisions)
sat_label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(
    data = sats_now, aes(x, y, label = label),
    family = font_mono, size = pt(5), colour = col_sub,
    segment.colour = NA, box.padding = 0.15, point.padding = 0.1,
    min.segment.length = Inf, seed = 1
  )
} else {
  txt(sats_now$x + inch(0.05), sats_now$y + inch(0.05), sats_now$label, 5)
}

p <- ggplot() +
  # --- sky ---------------------------------------------------
  geom_rect(data = horizon_glow,
            aes(xmin = cx0, xmax = cx1, ymin = ymin, ymax = ymax, alpha = alpha),
            fill = "#1E3A7A", colour = NA) +
  { if (draw_aurora)
      geom_rect(data = aurora,
                aes(xmin = xmin, xmax = xmax, ymin = y0, ymax = y1, fill = fill, alpha = alpha),
                colour = NA) } +
  geom_point(data = stars, aes(x, y, size = size, alpha = alpha),
             colour = col_star, stroke = 0) +
  geom_path(data = tracks, aes(x, y, group = grp),
            colour = col_orbit, linewidth = 0.18, linetype = "22") +
  geom_segment(data = signals, aes(x = x, y = y, xend = xend, yend = yend),
               colour = col_signal, linewidth = 0.12, alpha = 0.35) +
  # satellite glyph: body + solar panels
  geom_segment(data = sats_now, aes(x = x - inch(0.07), xend = x + inch(0.07), y = y, yend = y),
               colour = col_text, linewidth = 0.35) +
  geom_point(data = sats_now, aes(x, y), shape = 15, size = 1.1, colour = col_text) +
  geom_point(data = sats_now, aes(x, y), shape = 16, size = 4, colour = col_signal, alpha = 0.12) +
  sat_label_layer +
  annotate("segment", x = cx0 + safe, xend = cx1 - safe, y = horizon_y, yend = horizon_y,
           colour = col_edge, linewidth = 0.25) +
  txt(compass$x, horizon_y + inch(0.03), compass$lab, 6, hjust = 0.5) +
  # --- ground ------------------------------------------------
  geom_sf(data = graticule, colour = col_grat, linewidth = 0.15) +
  geom_sf(data = states, fill = col_land, colour = NA) +
  geom_sf(data = no_estimate, fill = NA, colour = col_edge, linewidth = 0.2) +
  pass_layers(passes) +
  geom_sf(data = headlands, aes(colour = colour), linewidth = 0.18, alpha = 0.45) +
  annotate("point", x = home[1], y = home[2], shape = 21, size = 1.6,
           fill = col_bg, colour = col_text, stroke = 0.4) +
  txt(home[1] + inch(0.06), home[2] - inch(0.02), "MHK", 5, colour = col_text, vjust = 1) +
  # --- title (sky, top left) --------------------------------
  txt(cx0 + safe, sky_top, "GUIDED FROM ABOVE", 26, family = font_title,
      colour = col_text, vjust = 1, fontface = "bold") +
  txt(cx0 + safe, sky_top - inch(0.45),
      "GPS on U.S. cropland, 2013\u20132019", 10, family = font_title, vjust = 1) +
  txt(cx1 - safe, sky_top,
      sprintf("GPS SKY OVER MANHATTAN, KS\n%s UTC%s",
              format(obs_time, "%Y-%m-%d %H:%M"),
              if (sky_is_real) "" else "\n(approximate)"),
      6, hjust = 1, vjust = 1) +
  # --- legend (bottom right) --------------------------------
  geom_sf(data = key_passes, colour = col_signal, linewidth = 0.16) +
  geom_sf(data = key_boxes, fill = NA, colour = col_edge, linewidth = 0.2) +
  txt(key_labels$x, key_y0 - inch(0.05), key_labels$lab, 6, hjust = 0.5, vjust = 1) +
  txt(lg_x0, key_y0 + key_side + inch(0.07),
      "PASS SPACING (denser = higher)\n2019 GPS share of cropland acres", 6.5, colour = col_text) +
  geom_rect(data = ramp_df, aes(xmin = xmin, xmax = xmax, ymin = ramp_y, ymax = ramp_y + ramp_h,
                                fill = fill), colour = NA) +
  # end labels are justified inward so they stay inside the safe area
  txt(ramp_ticks$x[2:4], ramp_y - inch(0.04), ramp_ticks$lab[2:4], 6, hjust = 0.5, vjust = 1) +
  txt(ramp_ticks$x[1],   ramp_y - inch(0.04), ramp_ticks$lab[1],   6, hjust = 0,   vjust = 1) +
  txt(ramp_ticks$x[5],   ramp_y - inch(0.04), ramp_ticks$lab[5],   6, hjust = 1,   vjust = 1) +
  txt(lg_x0, ramp_y + ramp_h + inch(0.07),
      "PASS COLOUR\nchange 2013\u20132019, pct. points", 6.5, colour = col_text) +
  txt(lg_x1, lg_y0, credit, 5.5, hjust = 1) +
  # --- scales & frame ---------------------------------------
  scale_colour_identity() +
  scale_fill_identity() +
  scale_alpha_identity() +
  scale_size_identity() +
  coord_sf(crs = st_crs(crs_conus), xlim = c(cx0, cx1), ylim = c(cy0, cy1),
           expand = FALSE, datum = NA) +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = col_bg, colour = NA),
    panel.background = element_rect(fill = col_bg, colour = NA),
    plot.margin      = margin(0, 0, 0, 0),
    legend.position  = "none"
  )

ggsave(out_file, p, width = canvas_w_in, height = canvas_h_in,
       dpi = export_dpi, bg = col_bg, device = "png")

# Optional proof: same image with the trim and safe lines drawn
if (use_showtext) showtext::showtext_opts(dpi = 600)   # must match the proof dpi
proof <- p +
  annotate("rect", xmin = cx0 + inch(0.25), xmax = cx1 - inch(0.25),
           ymin = cy0 + inch(0.25), ymax = cy1 - inch(0.25),
           fill = NA, colour = "red", linewidth = 0.3, linetype = "dashed") +
  annotate("rect", xmin = cx0 + safe, xmax = cx1 - safe,
           ymin = cy0 + safe, ymax = cy1 - safe,
           fill = NA, colour = "orange", linewidth = 0.3, linetype = "dotted")
ggsave("guided_from_above_PROOF.png", proof, width = canvas_w_in, height = canvas_h_in,
       dpi = 600, bg = col_bg, device = "png")
if (use_showtext) showtext::showtext_opts(dpi = export_dpi)

# ============================================================
# 8. WRITE-UP FOR THE SUBMISSION EMAIL (<= 75 words)
#    The Gannon claim is made only when the sky is built from
#    element sets near obs_time.
# ============================================================

sky_sentence <- if (sky_is_real) {
  paste("Above: the actual GPS constellation over Manhattan, Kansas during the May 2024",
        "Gannon geomagnetic storm, propagated from historical Space-Track orbital elements",
        "(two-body + J2) and unrolled into a panoramic sky.")
} else {
  paste("Above: an approximate GPS constellation over Manhattan, Kansas, propagated from",
        "orbital elements (two-body + J2) and unrolled into a panoramic sky.")
}

writeup <- paste(
  "Built entirely in R (sf, tigris, ggplot2). Each state is plowed with automated-guidance",
  "passes whose density reflects the 2019 share of cropland acres using GPS (USDA ERS);",
  "pass colour shows the 2013-2019 percentage-point change.",
  sky_sentence,
  "Census TIGER boundaries, CONUS Albers (EPSG:5070)."
)
n_words <- lengths(strsplit(writeup, "\\s+"))
cat(writeup, "\n\nWords:", n_words, "\n")
if (n_words > 75) warning("Write-up exceeds GeoHipster's 75-word limit.", call. = FALSE)
