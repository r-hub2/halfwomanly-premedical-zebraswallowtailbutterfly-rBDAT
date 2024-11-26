test_exception <- function() {
    RUnit::checkException(rBDAT:::throw("Hello, error!"))
}
