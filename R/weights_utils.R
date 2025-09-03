#' @title Configure a custom path for the neural-network weights (.rds)
#' @description
#' Sets (and validates) a custom file path for the pre-saved weights object
#' used by the package. This function does **not** read the file immediately;
#' it only stores the path and clears the in-memory cache so the next call to
#' [get_weights_dense()] will reload from disk.
#'
#' Package-default behavior looks for a file named
#' `inst/models/weights_nn_model.rds` at install time and resolves it via
#' [base::system.file()]. Use this function if you want to override that path.
#'
#' @param path `character(1)`. Absolute or relative path to a `.rds` file
#'   that exists and is readable.
#' @return Invisibly returns the path (invisible character scalar).
#' @examples
#' \dontrun{
#' # Point to a custom weights file:
#' set_weights_path("/path/to/weights_nn_model.rds")
#'
#' # Then use it (will be loaded lazily on first call):
#' W <- get_weights_dense()
#' }
#' @seealso [get_weights_dense()]
#' @export
set_weights_path <- function(path) {
  if (!file.exists(path)) stop("Path does not exist: ", path)
  assign("weights_path", path, envir = .astrolabe_cache)
  assign("weights_dense", NULL, envir = .astrolabe_cache)  # reset cache
  invisible(path)
}

#' @title Retrieve the neural-network weights (lazy-loaded, cached)
#' @description
#' Returns the weights object (previously saved with [base::saveRDS()])
#' used by functions in this package. The object is loaded on first call and
#' cached in memory for subsequent calls. If a custom path was set via
#' [set_weights_path()], that file is used. Otherwise, the function looks for
#' `inst/models/weights_nn_model.rds` bundled with the package using
#' [base::system.file()].
#'
#' @details
#' This function avoids reading files at top-level (package load time) to keep
#' `devtools::load_all()` and `devtools::document()` robust. The cache is
#' stored in a private environment and persists while the package is loaded.
#'
#' @return The R object read from the `.rds` file (typically a list/matrix
#'   of layer weights).
#' @examples
#' \dontrun{
#' # Default (uses the file inside the installed package, if present):
#' W <- get_weights_dense()
#'
#' # After setting a custom path:
#' set_weights_path("/path/to/weights_nn_model.rds")
#' W2 <- get_weights_dense()
#' }
#' @seealso [set_weights_path()]
#' @export
get_weights_dense <- function() {
  # Fast path: return from in-memory cache if available
  w <- get0("weights_dense", envir = .astrolabe_cache, inherits = FALSE)
  if (!is.null(w)) return(w)

  # Resolve path: user override via set_weights_path(), else inst/models/
  path <- get0("weights_path", envir = .astrolabe_cache, inherits = FALSE)
  if (is.null(path)) {
    path <- system.file("models", "weights_nn_model.rds", package = "astrolabe")
  }

  if (!nzchar(path) || !file.exists(path)) {
    stop("Could not find 'weights_nn_model.rds'. ",
         "Place it under inst/models/ in the package or call set_weights_path().")
  }

  # Read once, then cache for future calls (binary-safe across OSes)
  w <- readRDS(path)
  assign("weights_dense", w, envir = .astrolabe_cache)
  w
}

# ---- internal cache (not exported) ------------------------------------------

#' Internal package cache (environment)
#'
#' @keywords internal
#' @noRd
.astrolabe_cache <- new.env(parent = emptyenv())
