#' Calculate Genetic Gain Based on Different Selection Models
#'
#' This function estimates genetic gain using various selection models (e.g., HS, GBlup, ssGBlup, etc.) based on breeding values (BV), genomic values (GV), and heritability parameters.
#'
#' @param Dt A data frame containing the dataset to be analyzed. Default is `OP.SIMUL$`Self.0_nSNPs = 400_Rep1`$Data`.
#' @param BV. A data frame containing breeding values and phenotypic data. Default is `BV.OP`.
#' @param IS Selection intensity (proportion of top individuals to select). Default is `0.01`.
#' @param model A character vector specifying which models to use for gain calculation. Options include `"HS"`, `"GBlup"`, `"ssGBlup"`, `"GDBlup"`, `"ssGDBlup"`, `"ssGDBlup100"`, and `"Pheno"`.
#' @param use_ssGBLUP Logical. Whether to include single-step GBLUP models in the analysis. Default is `FALSE`.
#' @param Group_By A character vector specifying grouping variables for summarizing results. Default is `c("SelfRate", "nSNPs")`.
#' @param Herit A data frame containing heritability estimates for each model and type. Default is `OP.SIMUL$`Self.0_nSNPs = 400_Rep1`$Herit.All`.
#' @param RealParameters A data frame containing true parameter values (e.g., real heritabilities). Default is `OP.SIMUL$`Self.0_nSNPs = 400_Rep1`$RealParameters`.
#' @param BV.ssGBlup.model Optional data frame containing BV estimates from ssGBlup models. Default is `NULL`.
#'
#' @return A list containing:
#' \describe{
#'   \item{GainResults}{A data frame with realized and expected genetic gains for each model and gain type (BV or GV).}
#' }
#'
#' @details
#' The function computes genetic gain by selecting the top individuals based on BV or GV and calculating the expected and realized gain using both estimated and true heritabilities. It supports multiple models and integrates both additive and dominance effects.
#'
#' @examples
#' \dontrun{
#' results <- GenGainFunc(Dt = my_data, BV. = my_BV, model = c("HS", "GBlup"), IS = 0.05)
#' head(results$GainResults)
#' }
#'
#' @export
#
GenGainFunc <- function(Dt = OP.SIMUL$`Self.0_nSNPs = 400_Rep1`$Data,
                        BV. = BV.OP,
                        IS = 0.01,
                        model = c("HS", "GBlup"),
                        use_ssGBLUP = FALSE,
                        Group_By = c("SelfRate", "nSNPs"),
                        Herit = OP.SIMUL$`Self.0_nSNPs = 400_Rep1`$Herit.All,
                        RealParameters = OP.SIMUL$`Self.0_nSNPs = 400_Rep1`$RealParameters,
                        BV.ssGBlup.model = NULL) {
  # Helper function to calculate BV and GV gains
  calc_gain_combined <- function(data,
                                 col_BV,
                                 col_GV,
                                 h2a,
                                 h2g,
                                 real_h2a,
                                 real_h2g,
                                 model_name) {
    # Get the mean
    mu <- mean(BV.[["Pheno"]], na.rm = T)
    BV_gain <- data |>
      arrange(desc(!!sym(col_BV))) |>
      head(nrow(data) * IS) |>
      group_by(across(all_of(Group_By))) |>
      reframe(
        Re_Exp_Gain_BreedEQ = round(((
          mean(Pheno, na.rm = TRUE) - mu
        ) * h2a) / mu, 3),
        Re_Exp_Gain = round((mean(
          !!sym(col_BV), na.rm = TRUE
        ) / mu) - 1, 3),
        Realised_Gain_BreedEQ = round(((
          mean(Pheno, na.rm = TRUE) - mu
        ) * real_h2a) / mu, 3),
        Realised_Gain = round((mean(
          !!sym("u+bv"), na.rm = TRUE
        ) / mu) - 1, 3),
        TypeGain = "BV",
        Model = model_name,
        SI = IS
      )
    
    GV_gain <- data |>
      arrange(desc(!!sym(col_GV))) |>
      head(nrow(data) * IS) |>
      group_by(across(all_of(Group_By))) |>
      reframe(
        Re_Exp_Gain_BreedEQ = round(((
          mean(Pheno, na.rm = TRUE) - mu
        ) * h2g) / mu, 3),
        Re_Exp_Gain = round((mean(
          !!sym(col_GV), na.rm = TRUE
        ) / mu) - 1, 3),
        Realised_Gain_BreedEQ = round(((
          mean(Pheno, na.rm = TRUE) - mu
        ) * real_h2g) / mu, 3),
        Realised_Gain = round((mean(
          !!sym("u+gv"), na.rm = TRUE
        ) / mu) - 1, 3),
        TypeGain = "GV",
        Model = model_name,
        SI = IS
      )
    
    bind_rows(BV_gain, GV_gain)
  }
  
  # Extract mean phenotype
  mu <- mean(BV.[["Pheno"]], na.rm = TRUE)
  
  # Extract heritabilities
  get_h2 <- function(model_name, type) {
    Herit |> filter(TypePop == model_name &
                      h2.type == type) |> pull(Estimate)
  }
  
  # Prepare results list
  results_list <- list()
  
  for (mod in c(model)) {
    # True parameter
    h2a <- get_h2("True", "h2a")
    h2g <- get_h2("True", "h2g")
    real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
    real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
    Data.GenGain <- BV. |> mutate(`u+bv` = BV + mu, `u+gv` = BV + mu + Dom_eff)
    results_list[["True"]] <- calc_gain_combined(Data.GenGain,
                                                 "u+bv",
                                                 "u+gv",
                                                 h2a,
                                                 h2g,
                                                 real_h2a,
                                                 real_h2g,
                                                 "True")
    # Phenotypic value
    if (mod == "Pheno") {
      h2a <- get_h2("HS", "h2a")
      h2g <- get_h2("HS", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV. |> mutate(
        Pheno = Pheno,
        Pheno2 = Pheno,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(
        data = Data.GenGain,
        col_BV = "Pheno",
        col_GV = "Pheno2",
        h2a =  h2a,
        h2g = h2g,
        real_h2a = real_h2a,
        real_h2g = real_h2g,
        model_name = mod
      )
    }
    # Half-Sibs model
    if (mod == "HS") {
      h2a <- get_h2("HS", "h2a")
      h2g <- get_h2("HS", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV. |> mutate(
        `u+a` = BV.HS + mu,
        `u+g` = BV.HS + mu,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(Data.GenGain,
                                                "u+a",
                                                "u+g",
                                                h2a,
                                                h2g,
                                                real_h2a,
                                                real_h2g,
                                                mod)
    }
    # GBlup model
    if (mod == "GBlup") {
      h2a <- get_h2("GBlup", "h2a")
      h2g <- get_h2("GBlup", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV. |> mutate(
        `u+a.gblup` = BV.GBlup + mu,
        `u+g.gblup` = BV.GBlup + DomEff.GD + mu,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(Data.GenGain,
                                                "u+a.gblup",
                                                "u+g.gblup",
                                                h2a,
                                                h2g,
                                                real_h2a,
                                                real_h2g,
                                                mod)
    }
    # GDBlup model
    if (mod == "GDBlup") {
      h2a <- get_h2("GDBlup", "h2a")
      h2g <- get_h2("GDBlup", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV. |> mutate(
        `u+a.gdblup` = BV.GDBlup + mu,
        `u+g.gdblup` = BV.GDBlup + DomEff.GD + mu,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(Data.GenGain,
                                                "u+a.gdblup",
                                                "u+g.gdblup",
                                                h2a,
                                                h2g,
                                                real_h2a,
                                                real_h2g,
                                                mod)
    }
    # ssGBlup model
    if (mod == "ssGBlup" && use_ssGBLUP) {
      h2a <- get_h2("ssGBlup", "h2a")
      h2g <- get_h2("ssGBlup", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV.ssGBlup.model |> mutate(
        `u+ssGBlup.a` = BV.ssGBlup + mu,
        `u+ssGBlup.g` = BV.ssGBlup + mu + Dom.ssGBlup,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(Data.GenGain,
                                                "u+ssGBlup.a",
                                                "u+ssGBlup.g",
                                                h2a,
                                                h2g,
                                                real_h2a,
                                                real_h2g,
                                                mod)
    }
    # ssGDBlup model
    if (mod == "ssGDBlup" && use_ssGBLUP) {
      h2a <- get_h2("ssGDBlup", "h2a")
      h2g <- get_h2("ssGDBlup", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV.ssGBlup.model |> mutate(
        `u+ssGDBlup.a` = BV.ssGDBlup + mu,
        `u+ssGDBlup.g` = BV.ssGDBlup + mu + Dom.ssGDBlup,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(Data.GenGain,
                                                "u+ssGDBlup.a",
                                                "u+ssGDBlup.g",
                                                h2a,
                                                h2g,
                                                real_h2a,
                                                real_h2g,
                                                mod)
    }
    # ssGDBlup full model
    if (mod == "ssGDBlup100" && use_ssGBLUP) {
      h2a <- get_h2("ssGBlup100", "h2a")
      h2g <- get_h2("ssGBlup100", "h2g")
      real_h2a <- RealParameters |> filter(VarType == "h2a", Pop == "F1") |> pull(Estimates)
      real_h2g <- RealParameters |> filter(VarType == "h2g", Pop == "F1") |> pull(Estimates)
      Data.GenGain <- BV.ssGBlup.model |> mutate(
        `u+ssGDBlup100.a` = BV.HDBlup + mu,
        `u+ssGDBlup100.g` = BV.HDBlup + mu + DomEff.HD,
        `u+bv` = BV + mu,
        `u+gv` = BV + mu + Dom_eff
      )
      results_list[[mod]] <- calc_gain_combined(Data.GenGain,
                                                "u+ssGDBlup100.a",
                                                "u+ssGDBlup100.g",
                                                h2a,
                                                h2g,
                                                real_h2a,
                                                real_h2g,
                                                mod)
    }
  }
  
  # Combine results into a single table
  GainResults <- bind_rows(results_list)
  
  return(list(GainResults = GainResults))
}

