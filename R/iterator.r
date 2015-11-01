#' An IRanges iterator
#'
#' @param obj A \code{\linkS4class{IRanges}} object.
#' @param ... Additional arguments affecting the iterator.
#'
#' @return The iterator
#' @export
#' @importFrom iterators iter
#' @examples
#' iR <- iter(IRanges(1:3, 10:13))
#' nextElem(iR)
#' nextElem(iR)
#' nextElem(iR)
iter.IRanges <- function(obj, ...) {
  i <- 0L
  len <- length(obj)

  nextEl <- function() {
    if (i == len) {
      stop("StopIteration")
    }
    i <<- i + 1L
    obj[i]
  }

  hasNxt <- function() {
    if (i < len) {
      TRUE
    } else {
      FALSE
    }
  }

  structure(
    list(nextElem = nextEl, hasNext = hasNxt),
    class = c("iRanges", "abstractiter", "iter")
  )
}

#' Does this iterator have a next element
#'
#' @param obj An iterator object.
#' @param ... Additional arguments that are ignored.
#'
#' @return Logical value
#' @export
#' @examples
#' ##
hasNext <- function(obj, ...) {
  UseMethod("hasNext")
}

#' @export
hasNext.ihasNext <- function(obj, ...) {
  obj$hasNext()
}

#' @export
hasNext.iRanges <- function(obj, ...) {
  obj$hasNext()
}

#' Create an iterator with a hasNext method
#'
#' @param it An iterable object.
#'
#' @return An \code{ihasNext} iterator
#' @export
#' @examples
#' ###
ihasNext <- function(it) {
  if (!is.null(it$hasNext)) {
    return(it)
  }
  cache <- NULL
  has_next <- NA

  nextEl <- function() {
    if (!hasNxt()) {
      stop("StopIteration", call. = FALSE)
    }
    has_next <<- NA
    cache
  }

  hasNxt <- function() {
    if (!is.na(has_next)) {
      return(has_next)
    }
    tryCatch({
      cache <<- nextElem(it)
      has_next <<- TRUE
    },
    error = function(e) {
      if (identical(conditionMessage(e), "StopIteration")) {
        has_next <<- FALSE
      } else {
        stop(e)
      }
    })
    has_next
  }

  structure(
    list(nextElem = nextEl, hasNext = hasNxt),
    class = c("ihasNext", "abstractiter", "iter")
  )
}


