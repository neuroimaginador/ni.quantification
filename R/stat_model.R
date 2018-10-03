#' Statistical modelling for Connectome Parameters
#'
#' @param subject_demographic (age and gender)
#' @param subject_parameters contain result of computing connectome parameters for one subject (nodeDegree, centrality, pairsConnectivity,
#' betweenness, smallworldness, mean degree)
#'
#' @return
#'
stat_model <- function(subject_demographic, subject_parameters){

  # Load from database
  load(system.file("connectome_database.Rdata", package = "debra.algorithms.connectome"))
  # The database will contain two variables, demographic and parameters
  # demographic will be a data frame with columns age and gender
  # parameters is a list with as many components as rows in the demographic variable
  # and each component is a list of the same structure as the one in subject_parameters.
  # Thus, parameters are going to match by name.

  age <- demographic$age
  gender <- demographic$gender

  #filter database corresponding to the subject
  min_age <- subject_demographic$age - 5
  max_age <- subject_demographic$age + 5

  #Filter gender
  filter_gender <- which(gender == subject_demographic$gender)
  age_gender <- age[filter_gender]

  minimum_sample_size <- 20

  if (length(filter_gender) >= minimum_sample_size) {

    #Age interval (minimum 20 subjects in sample)
    k <- 0
    repeat {

      f <- subset(age_gender, age_gender > min_age - k & age_gender < max_age + k)
      k <- k+1

      if (length(f) >= 20) {break}
    }

    #Gender filtered database
    filtered_database_gender <- parameters[filter_gender]

    #Age filtered database
    filter_age <- which(age_gender <= max(f) & age_gender >= min(f))
    filtered_database <- filtered_database_gender[filter_age]

  } else {

    filtered_database <- parameters[filter_gender]

  }


  # Loop over all parameters, by name:
  param_names <- names(subject_parameters)

  result <- list()

  for (name in param_names) {

    subject_value <- subject_parameters[[name]]
    database_values <- lapply(filtered_database, function(s) s[[name]])
    n_params <- unique(unlist(lapply(database_values, function(s) length(s))))

    if (length(subject_value) == 1) {

      # Global parameter
      database_values <- as.vector(unlist(database_values))
      database_values <- database_values[!is.na(database_values)]

      if ((length(unique(database_values)) > 1) && (!is.na(subject_value))) {

        distribution_function <- ecdf(database_values)

        quantile_subject <- distribution_function(subject_value)

        outlier <- "none"
        if (quantile_subject < 0.05) {

          outlier <- "below"

        } else {

          if (quantile_subject > 0.95) {

            outlier <- "above"

          }

        }

      } else {

        outlier <- "none"

      }
      result[[name]] <- list(list(value = subject_value, outlier = outlier))

    } else {

      # Convert database values to matrix
      values <- t(matrix(unlist(database_values), nrow = n_params))

      res <- list()
      for (i in 1:n_params) {

        vector_vals <- values[, i]
        vector_vals <- vector_vals[!is.na(vector_vals)]

        if ((length(unique(vector_vals)) > 1) && (!is.na(subject_value[i]))) {

          distribution_function <- ecdf(vector_vals)

          quantile_subject <- distribution_function(subject_value[i])

          outlier <- "none"
          if (quantile_subject < 0.05) {

            outlier <- "below"

          } else {

            if (quantile_subject > 0.95) {

              outlier <- "above"

            }

          }

        } else {

          outlier <- "none"

        }

        res[[i]] <- list(value = subject_value[i], outlier = outlier)

      }

      result[[name]] <- res

    }


  }

  return(result)

}




