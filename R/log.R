## Recording of everything a montest() call emits: warnings, messages, errors and printed output.

## `entering` is TRUE between capture_log() re-issuing a call and that call starting.
.mt_log_state <- new.env(parent = emptyenv())

## Evaluates `cl` (a montest() call) in `env` and returns its value with the element `log`: a data.frame
## with one row per warning, message or chunk of printed output, in the order they were emitted. Nothing is
## muffled: everything still reaches the console and the caller's own handlers. If the call fails, the
## same data.frame, ending with the error, is attached to the error condition as `log`.
capture_log <- function(cl, env) {
  type <- character()
  text <- character()
  con <- rawConnection(raw(0L), "r+")
  n_read <- 0L
  sink(con, split = TRUE)
  on.exit({
    sink()
    close(con)
  })

  add <- function(ty, tx) {
    type[[length(type) + 1L]] <<- ty
    text[[length(text) + 1L]] <<- sub("\n+$", "", tx)
  }
  ## Output printed since the last call, so that it keeps its place among the conditions.
  add_output <- function() {
    bytes <- rawConnectionValue(con)
    if (length(bytes) > n_read) {
      add("output", rawToChar(bytes[(n_read + 1L):length(bytes)]))
      n_read <<- length(bytes)
    }
  }
  record <- function(ty, cnd) {
    add_output()
    add(ty, conditionMessage(cnd))
  }
  log_table <- function() data.frame(type = type, text = text, stringsAsFactors = FALSE)

  .mt_log_state$entering <- TRUE
  out <- withCallingHandlers(
    eval(cl, env),
    warning = function(w) record("warning", w),
    message = function(m) record("message", m),
    error = function(e) {
      record("error", e)
      e$log <- log_table()
      stop(e)
    }
  )
  add_output()
  out$log <- log_table()
  out
}
