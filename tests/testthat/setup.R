library(magrittr)

list.files(path = file.path("..", "..", "R", "helper-functions"), pattern = "\\.R$", full.names = TRUE) %>%
  lapply(source, echo = FALSE)

