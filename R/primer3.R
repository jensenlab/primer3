
#' @useDynLib primer3
#' @importFrom Rcpp sourceCpp
NULL

#' Thermodynamic calculations with Primer3
#' 
#' These functions allow direct access to the Primer3 thermodynamic libraries via Rcpp.
#' For all functions besides \code{calculate_tm}, Primer3 requires initialization to
#' load the thermodynamic configuration data. The initialization will happen automatically
#' (and only once for each session) when a function that needs initialization is called.
#' (Or, the user can initialize manually by calling \code{primer3_init}, although there is
#' no advantage to this.) To avoid memory leaks, users should call \code{primer3_free} at
#' the end of the program to free memory allocated during initialization. Calling 
#' \code{primer3_free} after every call to a function will slow performance, since the 
#' configuration files will need to be reloaded during subsequent calls.
#' 
#' @param oligos Character vector of oligos. Calling once with multiple oligos is faster 
#' than repeated calls with single oligos.
#' @param salt_conc Monovalent salt concentration in mM.
#' @param divalent_conc Divalent ion concentration in mM.
#' @param dntp_conc dNTP concentration in mM.
#' @param dna_conc DNA concentration in nM.
#' @param nn_max_len Tm for oligos up to this length will be calculated using a nearest-neighbor
#' method. Beyond this length, the Tm is extrapolated based on GC content.
#' @param tm_method Method for Tm calculations. Options include "Breslauer" and "SantaLucia". 
#' As with \code{salt_correction}, this package uses the recommended Primer3
#' method ("SantaLucia"), not the Primer3 default method ("Schildkraut").
#' @param salt_correction Method for salt correction calculations. Options include "Schildkraut", 
#' "SantaLucia", "Owczarzy". As with \code{tm_method}, this package uses the recommended Primer3
#' method ("SantaLucia"), not the Primer3 default method ("Schildkraut").
#' @param maxloop Maximum length of loop structures. Not available for \code{calculate_tm}.
#' @param temp_c Temperature in degrees Celsius. Not available for \code{calculate_tm}.
#' @param print_output If \code{TRUE}, print alignment or secondary structure to terminal. 
#' When \code{TRUE}, no output is returned. Not available for \code{calculate_tm}.
#' 
#' @return For \code{calculate_tm}, a numeric vector with the melting temperature [C] for each
#' input oligo. For the other functions, a named list of vectors indicating if a structure was 
#' found (\code{structure_found}); the changes in entropy (\code{ds}), enthalpy (\code{dh}), 
#' and Gibbs free energy (\code{dg}); and the alignment locations for each end (\code{align_end_1}
#' and \code{align_end_2}). Note that if \code{print_output} is \code{TRUE}, the functions return
#' \code{NULL}.
#' 
#' @seealso \code{\link{primer3_free}}
#' 
#' @name thermo
NULL

TM_METHODS <- c("Breslauer", "SantaLucia")
SALT_CORRECTION_METHODS <- c("Schildkraut", "SantaLucia", "Owczarzy")

#' @rdname thermo
#' @export
calculate_tm <- function(oligos, 
                     salt_conc=50.0, 
                     divalent_conc=0.0, 
                     dntp_conc=0.0, 
                     dna_conc=50.0,
                     nn_max_len=60,
                     tm_method="SantaLucia",
                     salt_correction="SantaLucia") {
  tm_method <- pmatch(tolower(tm_method), tolower(TM_METHODS)) - 1L
  salt_correction <- pmatch(tolower(salt_correction), tolower(SALT_CORRECTION_METHODS)) - 1L
  if (is.na(tm_method) || is.na(salt_correction)) {
    stop("Invalid Tm or salt correction method.")
  }
  call_seq_tm(oligos, salt_conc, divalent_conc, dntp_conc, dna_conc, as.integer(nn_max_len), tm_method, salt_correction)
}

#' Loading Primer3 configuration files
#' 
#' For all functions besides \code{calculate_tm}, Primer3 requires initialization to
#' load the thermodynamic configuration data. The initialization will happen automatically
#' (and only once for each session) when a function that needs initialization is called.
#' (Or, the user can initialize manually by calling \code{primer3_init}, although there is
#' no advantage to this.) To avoid memory leaks, users should call \code{primer3_free} at
#' the end of the program to free memory allocated during initialization. Calling 
#' \code{primer3_free} after every call to a function will slow performance, since the 
#' configuration files will need to be reloaded during subsequent calls.
#' 
#' @param config_path Path to Primer3 configuration files. These are installed by default with
#' this package. The path must not end with a "\code{/}".
#' 
#' @name init
NULL

#' @rdname init
#' @export
primer3_init <- function(config_path=system.file("extdata/primer3_config", package="primer3")) {
  invisible(call_thal_init(paste0(config_path, .Platform$file.sep)))
}

#' @rdname init
#' @export
primer3_free <- function() {
  call_thal_free()
}

thal <- function(oligo1, oligo2, 
                 alignment_type = 1L,
                 maxloop = 30L,
                 mv = 50.0,
                 dv = 0.0,
                 dntp = 0.0,
                 dna = 50.0,
                 temp_c = 37.0,
                 debug = FALSE,
                 temp_only = FALSE,
                 dimer = FALSE,
                 print_output = FALSE) {
  if (!is_thal_init()) {
    primer3_init()
  }
  temp <- temp_c + 273.15  # temp must be absolute
  print_output <- as.integer(print_output)
  results <- call_thal(oligo1, oligo2, as.integer(debug), alignment_type, maxloop, mv, dv, dntp, dna, temp, as.integer(temp_only), as.integer(dimer), as.integer(print_output))

  if (print_output) {
    # the return values are messed up when printing; call twice
    return(NULL)
  } else {
    return(results)
  }
}

#' @rdname thermo
#' @export
calculate_hairpin <- function(oligo, ...) {
  thal(oligo, oligo, ..., alignment_type = 4L)
}

#' @rdname thermo
#' @export
calculate_homodimer <- function(oligo, ...) {
  thal(oligo, oligo, ..., alignment_type = 1L)
}

#' @rdname thermo
#' @export
calculate_dimer <- function(oligo1, oligo2, ...) {
  thal(oligo1, oligo2, ..., alignment_type = 1L)
}

