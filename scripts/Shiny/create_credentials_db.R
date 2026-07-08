# Run this script ONCE to create (or recreate) the credentials database.
# Usage: Rscript scripts/Shiny/create_credentials_db.R
#
# IMPORTANT: fill in real passwords below before running — do NOT commit passwords here.
# The resulting credentials.sqlite file is gitignored and never committed.

library(shinymanager)

credentials = data.frame(
  user     = c("admin",    "guest",    "sebastien"),
  password = c("CHANGE_ME", "CHANGE_ME", "CHANGE_ME"),  # replace before running
  admin    = c(TRUE,        FALSE,       FALSE),
  comment  = c("Admin user", "Read-only", "Sebastien"),
  stringsAsFactors = FALSE
)

if (any(credentials$password == "CHANGE_ME"))
  stop("Set real passwords in this script before running it. Do not commit the passwords.")

db_path = file.path("data", "credentials.sqlite")
#dir.create(dirname(db_path), showWarnings = FALSE)

create_db(
  credentials_data = credentials,
  sqlite_path      = db_path,
  passphrase       = NULL   # set a passphrase string here to encrypt the DB
)

#
cat("Credentials database created at:", db_path, "\n")
cat("Users:", paste(credentials$user, collapse = ", "), "\n")
cat("Run this script again any time you need to reset or add users.\n")




