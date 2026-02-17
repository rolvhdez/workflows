# Helper functions -----
check_out <- function(x) if (!dir.exists(x)) dir.create(x)
fancy_process <- function(
  process,
  spinner_type = "simpleDotsScrolling",
  message = "Processing",
  ...
) {
  #' Creates an environment to be executed with a 
  #' spinner function to show progress for a process
  #'
  #' @param process
  #' @param message Message to be displayed when process is executed
  # Define the wrapper function
  wrapper <- function() {
    tryCatch({
      # Start message
      cli_process_start(message)
      # Execute the process
      result <- do.call(process, list(...))
      cli_process_done()
      return(result)
    }, error = function(e) {
      # Finish with error message
      cli_alert_danger(paste("Error:", e$message))
      stop(e)
    })
  }
  wrapper()
}
read_sumstats_file <- function(sumstats_path, chunk_size = 1000000) {
  #' Read the summary statistics file
  #' 
  #' By default read in batches of 1 million lines
  #' to reduce memory usage.
  #'
  #' @param sumstats_path
  #' @param chunk_size
  #'
  con <- file(sumstats_path, "r")
  df <- data.table::fread(text = readLines(con, n = chunk_size))
  while (TRUE) {
    chunk <- readLines(con, n = chunk_size)

    # When the number of lines left is zero, break
    if (length(chunk) == 0) break

    c <- data.table::fread(text = chunk)
    if (!identical(names(c), names(df))) {
      colnames(c) <- names(df)
    }
    df <- data.table::rbindlist(list(df, c))
  }
  close(con)
  return(df) # Return
}

reformat_sumstats <- function(sumstats, model) {
  #' Change the table format to follow the template
  #' from https://r-graph-gallery.com/101_Manhattan_plot.html
  # 
  #' @param sumstats Raw GWAS summary statistics file
  #' @param model Software where it comes from
  
  require(dplyr)
  require(stats)
   
  # Column mappings for all available models
  models <- list(
    snipar = list(
      filter_col = "direct_log10_P",
      select_cols = c("chromosome", "pos", "SNP", "direct_log10_P", "direct_Beta", "freq", "direct_N"),
      new_names = c("CHR", "BP", "SNP", "P", "BETA", "MAF", "N")
    ),
    regenie = list(
      filter_col = "LOG10P",
      select_cols = c("CHROM", "GENPOS", "ID", "P", "BETA", "N"),
      new_names = c("CHR", "BP", "SNP", "P", "BETA", "N")
    )
  )
  if (!model %in% names(models)) {
    stop("Supported models are: ", paste(names(models), collapse = ", "))
  }
  spec <- models[[model]]
  x <- sumstats %>%
    dplyr::filter(!is.na(.data[[spec$filter_col]])) %>%
    dplyr::select(dplyr::all_of(spec$select_cols)) %>%
    stats::setNames(spec$new_names)
  # Apply model-specific transformations
  if (model == "snipar") {
    x <- x %>% 
      dplyr::mutate(
        BETA = BETA * -1,
        P = 10^(-P),
        MAF = ifelse(MAF > 0.5, 1 - MAF, MAF)
      )
  }
  
  return(x)
}
export_plot <- function(plot_obj, file_path, 
                       width = 1080, height = 1080 * 0.75, 
                       res = 150, units = "px") {
  #' Export a plot to file
  #'
  #' @param plot_obj The plot object to export
  #' @param plot_name Name of the plot (without extension)
  #' @param output_dir Output directory path
  #' @param width Plot width
  #' @param height Plot height
  #' @param res Plot resolution
  #' @param units Plot units
  
  # Export the plot
  png(
    filename = file_path,
    width = width, height = height,
    res = res, units = units
  )
  suppressMessages(print(plot_obj))
  dev.off()
  
  # Success message
  cli::cli_alert_success("Exported `{file_path}`")
  
  # Return the file path invisibly
  invisible(file_path)
}

# Custom font ---
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(showtext))
font_add_google("Source Sans 3", "source-sans-3")

# Custom colors ---
# https://coolors.co/palette/ff595e-ffca3a-8ac926-1982c4-6a4c93
red <- "#ff595e"
yellow <- "#FFCA3A"
green <- "#8AC926"
blue <- "#1982c4"
purple <- "#6A4C93"

# Plot theme
theme_set(
  theme_bw() +
    theme(
      text = element_text(family = "Source Sans 3"),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black",
                                 linewidth = 0.5),
      axis.line.y = element_line(color = "black",
                                 linewidth = 0.5),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(color = "#3d3d3d", size = 8),
      plot.caption = element_text(color = "#3d3d3d", size = 8),
      strip.text = element_text(color = "#3d3d3d", face = "bold", size = 10),
      strip.background = element_rect(
        color = "#3d3d3d", fill = "white", linewidth = 1
      )
    )
)
