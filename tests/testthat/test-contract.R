test_that("contract lists required fields", {
  ct <- amr_contract()
  expect_true(all(c("isolate_id", "species", "drug") %in% ct$field[ct$requirement == "required"]))
  expect_false(anyDuplicated(ct$field) > 0)
})

test_that("JSON schema exists, parses and matches the contract", {
  skip_if_not_installed("jsonlite")
  schema <- jsonlite::fromJSON(amr_contract_schema(), simplifyVector = FALSE)
  expect_setequal(names(schema$properties), amr_contract()$field)
})
