#' canavg
#' This function serves two purposes:
#' 1) Computes points in the sky associated with all survey points on the surface of a concave hemispherical mirror densiometer
#'
#' @param gla_detrep raw gla detailed report output


canavg <- function(gla_detrep = NA, eyey=-8, eyez=10, survmeth = "card", weightscalc = TRUE,  diagnostics=TRUE) {

  ### load in packages
  packages <- c("car","ggplot2","rgl","scatterplot3d","tidyverse")
  lapply(packages, library, character.only = TRUE)

  ### Define constants and XY grid of points on mirror surface
  rad_to_deg_coeff <- 180 / pi    # Coefficient to convert radians to degrees
  deg_to_rad_coeff <- pi / 180    # Coefficient to convert degrees to radians
  r_sphere <- 3                   # Radius of sphere (inches), provided by mirror densiometer user manual
  c_sphere <- r_sphere * 2 * pi   # Circumference of sphere (inches)
  quarter_c_sphere <- c_sphere/4  # One quarter of the circumference (inches)

  if (survmeth == "card") {
    # name cardinal directions
    dirscomplete <- c("n","e","s","w")
    # X and Y coordinates of 96 engraved surface points, treating mirror surface as a flat plane w/o hemispherical shape
    eng_coords <- tibble::tibble(
      x_eng_vect = c(0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.0625,0.1875,0.3125,0.4375,0.0625,0.1875,0.3125,0.4375,0.0625,0.1875,0.0625,0.1875,-0.0625,-0.1875,-0.3125,-0.4375,-0.5625,-0.6875,-0.0625,-0.1875,-0.3125,-0.4375,-0.5625,-0.6875,-0.0625,-0.1875,-0.3125,-0.4375,-0.0625,-0.1875,-0.3125,-0.4375,-0.0625,-0.1875,-0.0625,-0.1875,0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.0625,0.1875,0.3125,0.4375,0.0625,0.1875,0.3125,0.4375,0.0625,0.1875,0.0625,0.1875,-0.0625,-0.1875,-0.3125,-0.4375,-0.5625,-0.6875,-0.0625,-0.1875,-0.3125,-0.4375,-0.5625,-0.6875,-0.0625,-0.1875,-0.3125,-0.4375,-0.0625,-0.1875,-0.3125,-0.4375,-0.0625,-0.1875,-0.0625,-0.1875),
      y_eng_vect = c(0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.3125,0.3125,0.3125,0.3125,0.4375,0.4375,0.4375,0.4375,0.5625,0.5625,0.6875,0.6875,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.3125,0.3125,0.3125,0.3125,0.4375,0.4375,0.4375,0.4375,0.5625,0.5625,0.6875,0.6875,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.1875,-0.1875,-0.1875,-0.1875,-0.1875,-0.1875,-0.3125,-0.3125,-0.3125,-0.3125,-0.4375,-0.4375,-0.4375,-0.4375,-0.5625,-0.5625,-0.6875,-0.6875,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.1875,-0.1875,-0.1875,-0.1875,-0.1875,-0.1875,-0.3125,-0.3125,-0.3125,-0.3125,-0.4375,-0.4375,-0.4375,-0.4375,-0.5625,-0.5625,-0.6875,-0.6875))
  } else if (survmeth == "rose") {
    dirscomplete <- c("n", "ne", "e", "se", "s", "sw", "w", "nw")
    x_eng_vect <- NA
    y_eng_vect <- NA
  } else if (survmeth != "card" & survmeth != "rose") {
    stop("Invalid value for argument 'survmeth'. Use 'card or rose' only.")
  }

  ### Initialize output data frames
  # Tibble for XYZ coordinates of each of 96 engraved points
  # This tibble empty now - will be populated row-wise in for loop below
  out <- tibble(
    x_orgPlane = numeric(0),
    y_orgPlane = numeric(0),
    z_orgPlane = numeric(0),
    # point class grouping factor
    group = factor("mirror")
  )

  # Define 3-dimensional bounding box for later visualization
  # note - these values are not important for analysis, only setting bounds of 3D plot
  bound_box <- tibble(
    x_orgPlane = c(-3, -3, 3, 3, -3, -3, 3, 3),
    y_orgPlane = c(-3, 3, -3, 3, -3, 3, -3, 3),
    z_orgPlane = c(-3, -3, -3, -3, 3, 3, 3, 3),
    # point class grouping factor - this used solely for visualization
    group = factor("boundingbox")
  )

  # Define geometry for origin of 3D Cartesian space
  sphere_center_3d <- tibble(
    x_orgPlane = c(0),
    y_orgPlane = c(0),
    z_orgPlane = c(0),
    # point class grouping factor - this used solely for visualization
    group = factor("center")
  )

  ### compute XYZ coordinates n 3D space for each of 96 points on the mirror surface
  # note - anything called 'eng' or 'engraved' refers to the mirror surface w/ its engraved grid
  # note - this is not necessary to re-compute with each call to this function, but it's included in order to be transparent about the conversion process between the engraved grid and 3-dimensional space
  for (i in seq_len(nrow(eng_coords))) {
    # pick row corresponding to
    x_engraved <- eng_coords[i,1]
    y_engraved <- eng_coords[i,2]

    # Note - all XYZ coordinates refer to an origin at the center of the sphere partially described by the mirror surface

    # Step 1: Calculate X component of point on XY plane in XYZ coordinate system
    # calculate angle (radians) below XY plane of line w/ end points at origin and point's X value on engraved grid's X axis
    # here "negalt" indicates "negative altitude", or radians below horizon
    negalt_x <- ((90 - (x_engraved / (c_sphere / 4) * 90)) * -1) * deg_to_rad_coeff

    # Use known angle and hypotenuse length (3 in, radius of sphere) to calculate the X component of the point in XYZ space
    x_cartPlane <- cos(negalt_x) * r_sphere

    # Step 2: Calculate Y component of point on XY plane in XYZ coordinate system
    # calculate angle (radians) below XY plane of line w/ end points at origin and point's Y value on engraved grid's X axis
    # here "negalt" indicates "negative altitude", or radians below horizon
    negalt_y <- ((90 - (y_engraved / (c_sphere / 4) * 90)) * -1) * deg_to_rad_coeff
    # Use known angle and hypotenuse length (3 in, radius of sphere) to calculate the Y component of the point in XYZ space
    y_cartPlane <- cos(negalt_y) * r_sphere

    # Step 3: Calculate Z component of point in XYZ coordinate system
    # Calculate angle (radians) between line extending from origin to positive infinity on the XYZ Y axis and line with end points at origin and XY coordinates in XYZ space as calculated above
    az_ang_rad <- atan(x_cartPlane / y_cartPlane)
    # Calculate z value of sphere surface at XY point in XYZ space
    hyp <- sqrt(x_cartPlane^2 + y_cartPlane^2)
    z_cartPlane <- sqrt(r_sphere^2 - hyp^2) * -1

    # Store calculated values
    out[i,1] <- x_cartPlane
    out[i,2] <- y_cartPlane
    out[i,3] <- z_cartPlane
    # include string "mirror" to aid in later categorization and visualization
    out[i,4] <- "mirror"

    # When last iteration is ending, combine tibbles containing computed XYZ coordinates and bounding box and center coordinates
    if (i == nrow(eng_coords)) {
      combined_data <- na.omit(out) %>%
        bind_rows(bound_box) %>%
        bind_rows(sphere_center_3d)
    }
  }

  ### Compute slope and bearing of line reflecting off mirror surface emitted from eye
  # define point of light emission in XYZ space (eye position, with Y and X components defined in function call
  # note - these variables are given default values in function skypoints, but can be modified manually
  eye_position_xyz <- c(0, eyey, eyez)

  # Calculate plane tangent to sphere surface at each of 96 points
  # loop j is for stepping through each measurement direction
  for (j in 1:length(dirscomplete)) {
    if (j == 1) {
      # define empty list for storage of outputs
      list_out <- vector(mode='list', length=4)
      # code below defaults to E as azimuth 0, this corrects to appropriate cardinal direction correction such that due
      if (survmeth == "card") {
        correction <- c(-90,0,90,180)
      }
    }

    # note - i loop is nested within j loop
    # note - object out contains XYZ coordinates of 96 points
    for (i in 1:nrow(out)) {
      if (i == 1) {
        sphere_origin <- c(0, 0, 0)
        out2 <- out %>%
          mutate(azmAng = NA) %>%
          mutate(altAng = NA) %>%
          mutate(index = as.factor(dirscomplete[j])) %>%
          select(index, everything())
      }

      # Define the points
      reflection_point <- out2[i,2:4] %>% unlist() %>% as.vector()

      # Calculate the incident vector
      incident_vector <- reflection_point - eye_position_xyz

      # Calculate the normal vector at the reflection point
      # Since the mirror is a concave spherical mirror, the normal vector is just the normalized reflection point vector
      normal_vector <- reflection_point / sqrt(sum(reflection_point^2))

      # Calculate the reflection vector
      reflection_vector <- incident_vector - 2 * sum(incident_vector * normal_vector) * normal_vector

      # Calculate the azimuthal angle relative to the X-axis
      azimuthal_angle <- atan2(reflection_vector[2], reflection_vector[1]) * (180 / pi)
      if (azimuthal_angle < 0) {azimuthal_angle <- azimuthal_angle + 180}
      azimuthal_angle2 <- azimuthal_angle + correction[j]
      if (azimuthal_angle2 < 0) {azimuthal_angle2 <- azimuthal_angle2 + 360}
      out2[i,"azmAng"] <- azimuthal_angle2

      # Calculate the altitude above the horizon
      altitude <- asin(reflection_vector[3] / sqrt(sum(reflection_vector^2))) * (180 / pi)
      out2[i,"altAng"] <- altitude
    }
    list_out[[j]] <- out2
    if (j == 4) {
      out3 <- bind_rows(list_out) %>%
        mutate(total_brightness = NA_real_)
    }
  }

  # if using canavg to calculate weighted avg coefficients
  if (weightscalc == TRUE) {
  ### read in + use report from Gap Light Analyzer --------------------------------------------------
  # note - This is report output by GLA as .txt file, converted to .xlsx and
  # define column names by reading in .txt and pulling in needed strings
  raw_data <- readLines(gla_detrep)
  column_names_t1 <-  strsplit(raw_data[2], ";")[[1]] %>% gsub(pattern = " ", replacement = "", x = .)
  column_names_t2 <-  strsplit(raw_data[3], ";")[[1]] %>% gsub(pattern = " ", replacement = "", x = .)
  column_names <- paste0(column_names_t1,column_names_t2)
  rm(raw_data,column_names_t1,column_names_t2)

  # read in and manipulate main report
  gla_report <- read.delim(gla_detrep, sep = ";", header = FALSE, skip = 2) %>%
    # convert to tibble
    as_tibble() %>%
    # remove top row, which was relict of early data cleaning
    slice(-1) %>%
    # set column names
    set_names(column_names) %>%
    # convert columns 2-last from character to numeric
    mutate(across(.cols = 2:last_col(), .fns = ~parse_number(as.character(.)))) %>%
    # compute azimuth and altitude bin maxima degrees from bin number
    mutate(AziDegMax = AziBin*(360/length(unique(AziBin)))) %>%
    mutate(AltDegMax = AltBin*(90/length(unique(AltBin))))  %>%
    # compute azimuth and altitude bin minima degrees from bin number
    mutate(AziDegMin = (AziBin-1)*(360/length(unique(AziBin)))) %>%
    mutate(AltDegMin = (AltBin-1)*(90/length(unique(AltBin)))) %>%
    # compute azimuth and altitude bin means degrees from bin number
    mutate(AziDegCenter = (AziDegMax + AziDegMin)/2 ) %>%
    mutate(AltDegCenter = (AltDegMax + AltDegMin)/2 )

  # if incorrect number of azimuth and altitude bands used, throw error
  if (length(unique(gla_report$AziBin)) != 90 | length(unique(gla_report$AltBin)) != 45) {
    stop("Insufficient number of azimuth or altitude bins. \nEnsure that GLA output is divided into exactly 90 azimuth bins and 45 altitude bins.\nYou will need to return to GLA, reconfigure the output, and re-run calculations there.")}

  # apply gla_report$AboveTotal values to out3 rows based on azimuth and altitude
  # note - this is the important densiometer-GLA calibration step! All prior to this is preparation.
  for (k in 1:nrow(out3)) {
    # pull out individual azimuth and altitude values and numbers
    mirrorptAzi <- as.vector(unlist(out3[k,6]))
    mirrorptAlt <- as.vector(unlist(out3[k,7]))

    # pick azimuth-altitude bin associated w/ each point and extract insolation value
    out3[k,8] <- gla_report %>%
      # for azimuth - pick the column of azimuth bins which contains the point
      filter(AziDegMax > mirrorptAzi & AziDegMin < mirrorptAzi) %>%
      # for altitude - pick the row within this column which contains the point - this yields a single azi-alt bin
      filter(AltDegMax > mirrorptAlt & AltDegMin < mirrorptAlt) %>%
      # extract insolation value
      dplyr::select(AboveTotal)
  }

  #
  for (p in 1:4) {
    if (p == 1) {
      dirTotals <- rep(NA,4)
      names(dirTotals) <- dirscomplete
    }
    dirTotals[p] <- out3 %>%
      filter(index == dirscomplete[p]) %>%
      select(total_brightness) %>%
      sum()
  }

  ### print averaging weights and optionally print diagnostic plots ---------------------------------
  # compute weights
  dirWeights <- dirTotals/sum(dirTotals)

  # optionally display diagnostics
  if (diagnostics == TRUE) {
    # Generate 3D scatter plot for visual validation of geometry #
    car::scatter3d(
      y_orgPlane ~ x_orgPlane + z_orgPlane | group,
      data = combined_data,
      xlim = c(0, .2),
      ylim = c(0, 3),
      zlim = c(-3, 0),
      axis.ticks = TRUE,
      point.col = "black",
      sphere.size = 0.5,
      surface.alpha = 0,
      fogtype = "linear",
      residuals = FALSE
    )
    # Add main title
    title3d(main = "Mirror Surface Point in 3D Space", line = 2)
    # Add axis labels
    rgl::mtext3d(text = "X Axis (inches) - Perpindicular to Direction of View", edge = "x-", line = 3)
    rgl::mtext3d(text = "Y Axis (inches) - Parallel to Direction of View", edge = "y-", line = 3)
    rgl::mtext3d(text = "Z Axis (inches) - Normal Line to the Ground", edge = "z-", line = 3)

    # display sunlight contribution by each azimuth-altitude bin in arbitary unit
    dev.new()
    gg1 <- ggplot(gla_report, aes(x = AziDegCenter, y = AltDegCenter, fill = AboveTotal)) +
      ggtitle("Contribution of Sun by Azimuth-Altitude Bins") +
      geom_raster(alpha = 0.9) +  # Set raster transparency
      scale_fill_viridis_c() +  # Apply viridis color scale
      labs(
        x = "Azimuth\n(degrees clockwise from north)",
        y = "Altitude\n(degrees above horizon)",
        fill = "Legend\nUnits Determined\nin GLA",
        title = "Contribution of Sun by Azimuth-Altitude Bins"
      ) +
      scale_x_continuous(
        limits = c(0, 360),
        breaks = seq(0, 360, by = 10)
      ) +
      scale_y_continuous(
        limits = c(0, 90),
        breaks = seq(0, 90, by = 10)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "black", linewidth = 0.5),  # Black major gridlines
        panel.grid.minor = element_line(color = "black", linewidth = 0.25)  # Black minor gridlines
      )
    print(gg1)

    # visualize sky points corresponding to 96 mirror surface points
    dev.new()
    gg2<- ggplot(data = out3, aes(x = azmAng, y = altAng, color = factor(index))) +
      geom_point(shape = 19, size = 3) +  # Specify point shape
      scale_x_continuous(
        limits = c(0, 360),
        breaks = seq(0, 360, by = 10),
        name = "Azimuth\n(degrees clockwise from north)"  # X-axis label
      ) +
      scale_y_continuous(
        limits = c(0, 90),
        breaks = seq(0, 90, by = 10),
        name = "Altitude\n(degrees above horizon)"  # Y-axis label
      ) +
      labs(
        title = "Azimuth and altitude of each of 96 points in each of four cardinal directions",
        color = "Index"  # Legend title for color scale
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        plot.title = element_text(hjust = 0.5),  # Center plot title
        legend.position = "none",  # Remove legend
        panel.grid.major = element_line(color = "black", size = 0.5),  # Black major gridlines
        panel.grid.minor = element_line(color = "black", size = 0.25)  # Black minor gridlines
      )
    print(gg2)
    }
  # define output list
  # list contains all mirror XYZ coordinates, azimuth and altitude of assoc. points in sky, and insolation values

    out4 <- out3 %>% select(-group)

    output_full <- list(dirWeights,
                        out4,
                        gla_report,
                        "README: GLA detailed report passed to canavg, weights stored as output list [[1]]")

    names(output_full) <- c("direction_weights_vector","full_points_data_tibble","Raw GLA report","README")
  } else if (weightscalc == FALSE) {

    out4 <- out3 %>% select(-group)
    dirWeights <- NA
    gla_report <- NA
    output_full <- list(dirWeights,
                        out4,
                        gla_report,
                        "README: No GLA detailed report passed to canavg, weights not calculated")

    names(output_full) <- c("direction_weights_vector","full_points_data_tibble","Raw GLA report","README")
  }
  return(output_full)
}



