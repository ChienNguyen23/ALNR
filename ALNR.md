ĐỢT 2

- **Table**
    
    ```r
    # ==============================================================================
    # END-TO-END (ONE TABLE): Alb/NLR quartiles (Q1–Q4) as GROUPING
    # - ONE Word table only
    # - MASLD and SLF included as VARIABLES (rows)
    #   * MASLD: computed on analytic set (Alb/NLR available + masld_group defined)
    #   * SLF: computed ONLY among MASLD with valid slf_group (others NA; denominator varies by quartile)
    # - FIX: SLF is NOT dropped even if missingness > 50%
    # - Prints: exclusion counts + remaining N + Alb/NLR cutpoints + quartile sizes
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(officer)
        library(flextable)
        library(rlang)
    })
    
    # ----------------------------- 0) HELPERS -----------------------------
    ensure_cols <- function(df, cols, default = NA) {
        miss <- setdiff(cols, names(df))
        if (length(miss)) for (m in miss) df[[m]] <- default
        df
    }
    
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    fix_border_issues <- function(ft) {
        ft <- border_remove(ft)
        ft <- hline_top(ft, border = fp_border(color = "black", width = 1), part = "all")
        ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "all")
        ft <- hline(ft, border = fp_border(color = "grey70", width = 0.5), part = "body")
        ft
    }
    
    fmt_cut <- function(x, digits = 3) format(round(x, digits), nsmall = digits)
    
    # ----------------------------- 1) ALCOHOL -----------------------------
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_cols(data, c("ALQ121", "ALQ130"), default = NA)
        data %>%
            mutate(
                ALQ121 = safe_num(ALQ121),
                ALQ130 = safe_num(ALQ130),
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # ----------------------------- 2) EXCLUSION (WITH LOG) -----------------------------
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ ---\n")
        initial_rows_total <- nrow(data)
        
        data <- ensure_cols(
            data,
            c("RIAGENDR","RIDAGEYR","LUXSMED","LUXCAPM",
              "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO",
              "LBXSAL","LBDSALSI",
              "ALQ121","ALQ130",
              "RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV"),
            default = NA
        )
        
        log_tbl <- tibble::tibble(step = character(), before = integer(), after = integer(), removed = integer())
        add_log <- function(step, before, after) {
            log_tbl <<- bind_rows(log_tbl, tibble::tibble(
                step = step, before = before, after = after, removed = before - after
            ))
        }
        
        # 1) Viral hepatitis
        cat("1. Loại bỏ viêm gan virus...\n")
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR", "LBDHCV")
        hep_vars_exist <- intersect(hep_vars, names(data))
        before_n <- nrow(data)
        if (length(hep_vars_exist) > 0) {
            data <- data %>%
                mutate(across(all_of(hep_vars_exist), safe_num)) %>%
                filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
        }
        after_n <- nrow(data)
        add_log("Viral hepatitis", before_n, after_n)
        cat("-> Đã loại bỏ", before_n - after_n, "hàng.\n\n")
        
        # 2) Main measurements
        cat("2. Loại bỏ thiếu FibroScan/Albumin/Neut-Lymph...\n")
        before_n <- nrow(data)
        
        data <- data %>% drop_na(LUXSMED, LUXCAPM)
        
        data <- data %>%
            mutate(
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                Albumin_g_dL_tmp = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                )
            ) %>%
            filter(!is.na(Albumin_g_dL_tmp)) %>%
            mutate(
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXLYPCT = safe_num(LBXLYPCT),
                has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO),
                has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT)
            ) %>%
            filter(has_counts | has_pct)
        
        after_n <- nrow(data)
        add_log("Missing FibroScan/Albumin/Neut-Lymph", before_n, after_n)
        cat("-> Đã loại bỏ", before_n - after_n, "hàng.\n\n")
        
        # 3) Alcohol threshold
        cat("3. Loại bỏ rượu bia vượt ngưỡng...\n")
        data <- calculate_alcohol_consumption(data)
        before_n <- nrow(data)
        data <- data %>%
            mutate(RIAGENDR = safe_num(RIAGENDR)) %>%
            filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        after_n <- nrow(data)
        add_log("Alcohol above threshold", before_n, after_n)
        cat("-> Đã loại bỏ", before_n - after_n, "hàng.\n\n")
        
        # 4) Age >= 18
        cat("4. Loại bỏ tuổi < 18...\n")
        before_n <- nrow(data)
        data <- data %>%
            mutate(RIDAGEYR = safe_num(RIDAGEYR)) %>%
            filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        after_n <- nrow(data)
        add_log("Age < 18", before_n, after_n)
        cat("-> Đã loại bỏ", before_n - after_n, "hàng.\n\n")
        
        cat("--- HOÀN TẤT LOẠI TRỪ ---\n")
        cat("Tổng ban đầu:", initial_rows_total, "\n")
        cat("Còn lại sau loại trừ:", nrow(data), "\n")
        cat("Tổng đã loại:", initial_rows_total - nrow(data), "\n\n")
        
        list(data = data, log = log_tbl)
    }
    
    # ----------------------------- 3) DIAGNOSTIC (MASLD/SLF) -----------------------------
    define_diagnostic_criteria <- function(data) {
        cat("Bắt đầu: Xác định MASLD/SLF\n")
        
        data <- ensure_cols(
            data,
            c("LUXCAPM","LUXSMED","RIAGENDR","BMXBMI","BMXWAIST",
              "LBXGLU","LBXGH","DIQ010","DIQ050","DIQ070",
              "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
              "BPQ040A","LBDTRSI","BPQ090D","LBDHDDSI",
              "weekly_alcohol_grams"),
            default = NA
        )
        
        out <- data %>%
            mutate(
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, safe_num(LUXCAPM) >= 263),
                
                risk_bmi_waist = if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)), NA,
                    (safe_num(BMXBMI) >= 25) |
                        (safe_num(RIAGENDR) == 1 & safe_num(BMXWAIST) >= 94) |
                        (safe_num(RIAGENDR) == 2 & safe_num(BMXWAIST) >= 80)
                ),
                
                risk_glucose_diabetes = if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBXGLU), FALSE, safe_num(LBXGLU) >= 100)) |
                        (if_else(is.na(LBXGH), FALSE, safe_num(LBXGH) >= 5.7)) |
                        (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                
                any_bp_measurement_high =
                    (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, safe_num(BPXOSY1) >= 130 | safe_num(BPXODI1) >= 85)) |
                    (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, safe_num(BPXOSY2) >= 130 | safe_num(BPXODI2) >= 85)) |
                    (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, safe_num(BPXOSY3) >= 130 | safe_num(BPXODI3) >= 85)),
                
                risk_blood_pressure = if_else(
                    is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
                    NA,
                    any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                
                risk_triglycerides = if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBDTRSI), FALSE, safe_num(LBDTRSI) >= 1.70)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                risk_low_hdl = if_else(
                    is.na(RIAGENDR) | is.na(LBDHDDSI),
                    NA,
                    ((safe_num(RIAGENDR) == 1 & safe_num(LBDHDDSI) < 1.0) |
                         (safe_num(RIAGENDR) == 2 & safe_num(LBDHDDSI) < 1.3)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                )
            ) %>%
            rowwise() %>%
            mutate(
                num_cardiometabolic_risks = sum(
                    c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
                    na.rm = TRUE
                ),
                has_one_plus_cardiometabolic_risk = if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA,
                    num_cardiometabolic_risks > 0
                )
            ) %>%
            ungroup() %>%
            mutate(
                is_light_drinker_final = case_when(
                    is.na(weekly_alcohol_grams) ~ NA,
                    safe_num(RIAGENDR) == 1 & weekly_alcohol_grams < 210 ~ TRUE,
                    safe_num(RIAGENDR) == 2 & weekly_alcohol_grams < 140 ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                
                masld_group = case_when(
                    has_hepatic_steatosis == FALSE ~ "non-MASLD",
                    has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ case_when(
                        is.na(is_light_drinker_final) ~ "MASLD",
                        is_light_drinker_final == TRUE ~ "MASLD",
                        is_light_drinker_final == FALSE ~ "non-MASLD",
                        TRUE ~ NA_character_
                    ),
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    TRUE ~ NA_character_
                ),
                
                slf_group = case_when(
                    masld_group != "MASLD" ~ NA_character_,
                    is.na(LUXSMED) ~ NA_character_,
                    safe_num(LUXSMED) >= 8.0 ~ "SLF",
                    safe_num(LUXSMED) < 8.0 ~ "non-SLF",
                    TRUE ~ NA_character_
                )
            )
        
        cat("Hoàn tất: MASLD/SLF\n\n")
        out
    }
    
    # ----------------------------- 4) BLOOD INDICES (NO WBC) -----------------------------
    define_blood_indices <- function(data) {
        cat("Bắt đầu: Tạo Albumin/Neutrophil/Lymphocyte, NLR, Alb/NLR\n")
        
        data <- ensure_cols(data, c("LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"), default = NA)
        
        out <- data %>%
            mutate(
                LBXLYPCT = safe_num(LBXLYPCT),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBDNENO  = safe_num(LBDNENO),
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                
                Albumin_g_dL = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                Albumin_g_L = case_when(
                    !is.na(LBDSALSI) ~ LBDSALSI,
                    is.na(LBDSALSI) & !is.na(LBXSAL) ~ LBXSAL * 10,
                    TRUE ~ NA_real_
                ),
                
                Neutrophil_percent = LBXNEPCT,
                Lymphocyte_percent = LBXLYPCT,
                Neutrophil_count   = LBDNENO,
                Lymphocyte_count   = LBDLYMNO,
                
                NLR = case_when(
                    !is.na(Neutrophil_count) & !is.na(Lymphocyte_count) & Lymphocyte_count > 0 ~
                        Neutrophil_count / Lymphocyte_count,
                    (is.na(Neutrophil_count) | is.na(Lymphocyte_count)) &
                        !is.na(Neutrophil_percent) & !is.na(Lymphocyte_percent) & Lymphocyte_percent > 0 ~
                        Neutrophil_percent / Lymphocyte_percent,
                    TRUE ~ NA_real_
                ),
                
                Alb_NLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0, Albumin_g_dL / NLR, NA_real_)
            )
        
        cat("Hoàn tất: Blood indices\n\n")
        out
    }
    
    # ----------------------------- 5) CATEGORICAL VARIABLES -----------------------------
    define_categorical_variables <- function(data) {
        data <- ensure_cols(
            data,
            c("LBXSCR","LBDSCRSI","RIDAGEYR","RIAGENDR","RIDRETH1",
              "BMXHIP","BMXWAIST",
              "DMDEDUC2","DMDMARTZ","SMQ020","SMQ040",
              "DIQ010","DIQ050","DIQ070","LBXGH","LBXSGL",
              "BPQ020","BPQ040A",
              "KIQ022","KIQ025",
              "MCQ220","MCQ160B","MCQ160C","MCQ160D","MCQ160E"),
            default = NA
        )
        
        data %>%
            mutate(
                creatinine_mg_dl = case_when(
                    "LBXSCR" %in% names(.) ~ safe_num(LBXSCR),
                    "LBDSCRSI" %in% names(.) ~ safe_num(LBDSCRSI) / 88.4,
                    TRUE ~ NA_real_
                ),
                RIDAGEYR = safe_num(RIDAGEYR),
                RIAGENDR = safe_num(RIAGENDR),
                RIDRETH1 = safe_num(RIDRETH1),
                BMXHIP   = safe_num(BMXHIP),
                BMXWAIST = safe_num(BMXWAIST),
                
                eGFR = if_else(
                    !is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                    175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                        if_else(RIAGENDR == 2, 0.742, 1) * if_else(RIDRETH1 == 4, 1.212, 1),
                    NA_real_
                ),
                Waist_Hip_Ratio = if_else(BMXHIP > 0, BMXWAIST / BMXHIP, NA_real_),
                
                Gender = factor(case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female"),
                                levels = c("Male", "Female")),
                
                Education_Level = factor(case_when(
                    DMDEDUC2 == 1 ~ "Less than 9th grade",
                    DMDEDUC2 == 2 ~ "9-11th grade",
                    DMDEDUC2 == 3 ~ "High school graduate/GED",
                    DMDEDUC2 == 4 ~ "Some college or AA degree",
                    DMDEDUC2 == 5 ~ "College graduate or above",
                    TRUE ~ NA_character_
                ), levels = c("Less than 9th grade", "9-11th grade", "High school graduate/GED",
                              "Some college or AA degree", "College graduate or above")),
                
                Marital_Status = factor(case_when(
                    DMDMARTZ == 1 ~ "Married",
                    DMDMARTZ == 2 ~ "Widowed",
                    DMDMARTZ == 3 ~ "Divorced",
                    DMDMARTZ == 4 ~ "Separated",
                    DMDMARTZ == 5 ~ "Never married",
                    DMDMARTZ == 6 ~ "Living with partner",
                    TRUE ~ NA_character_
                ), levels = c("Married", "Living with partner", "Never married", "Divorced", "Widowed", "Separated")),
                
                Race_Ethnicity = factor(case_when(
                    RIDRETH1 == 1 ~ "Mexican American",
                    RIDRETH1 == 2 ~ "Other Hispanic",
                    RIDRETH1 == 3 ~ "Non-Hispanic White",
                    RIDRETH1 == 4 ~ "Non-Hispanic Black",
                    RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
                    TRUE ~ NA_character_
                )),
                
                Smoking_Status = factor(case_when(
                    SMQ020 == 2 ~ "Non-smoker",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
                    TRUE ~ NA_character_
                ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
                
                Diabetes_Mellitus = factor(
                    if_else(
                        DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                            (!is.na(LBXGH) & safe_num(LBXGH) >= 6.5) |
                            (!is.na(LBXSGL) & safe_num(LBXSGL) >= 126),
                        "Yes", "No", missing = "No"
                    ),
                    levels = c("No", "Yes")
                ),
                
                Hypertension = factor(
                    if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"),
                    levels = c("No", "Yes")
                ),
                
                Chronic_Kidney_Disease = factor(
                    if_else((KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60)), "Yes", "No", missing = "No"),
                    levels = c("No", "Yes")
                ),
                
                Cancer = factor(if_else(MCQ220 == 1, "Yes", "No", missing = "No"),
                                levels = c("No", "Yes")),
                
                Cardiovascular_Disease = factor(
                    if_else(MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1, "Yes", "No", missing = "No"),
                    levels = c("No", "Yes")
                )
            )
    }
    
    # ----------------------------- 6) PIPELINE WRAPPER -----------------------------
    create_all_features <- function(data) {
        cat("--- BẮT ĐẦU: TẠO BIẾN ---\n")
        out <- data %>%
            define_diagnostic_criteria() %>%
            define_blood_indices() %>%
            define_categorical_variables()
        cat("--- HOÀN TẤT TẠO BIẾN ---\n\n")
        out
    }
    
    # ----------------------------- 7) MISSING REPORT -----------------------------
    generate_missing_report <- function(data) {
        data.frame(
            Variable = names(data),
            Missing_Count = sapply(data, function(x) sum(is.na(x))),
            Missing_Percentage = sapply(data, function(x) mean(is.na(x)))
        ) %>%
            filter(Missing_Count > 0) %>%
            mutate(Missing_Percentage = paste0(format(round(Missing_Percentage * 100, 2), nsmall = 2), "%")) %>%
            arrange(desc(Missing_Count)) %>%
            rename("Biến" = Variable, "Số lượng thiếu" = Missing_Count, "Tỷ lệ thiếu" = Missing_Percentage)
    }
    
    # ----------------------------- 8) QUARTILES (Alb/NLR) -----------------------------
    compute_albnlr_cuts <- function(df, var = "Alb_NLR", probs = c(0.25, 0.50, 0.75)) {
        v <- df[[var]]
        v <- v[is.finite(v) & !is.na(v)]
        if (length(v) < 4) return(c(NA_real_, NA_real_, NA_real_))
        as.numeric(quantile(v, probs = probs, na.rm = TRUE, type = 2))
    }
    
    add_albnlr_quartiles <- function(df, cuts, var = "Alb_NLR", out_var = "AlbNLR_Q", digits = 3) {
        a <- cuts[1]; b <- cuts[2]; c <- cuts[3]
        if (any(is.na(cuts))) {
            df[[out_var]] <- NA
            return(df)
        }
        
        lab1 <- paste0("Q1 \u2264 ", fmt_cut(a, digits))
        lab2 <- paste0("Q2 (", fmt_cut(a, digits), "–", fmt_cut(b, digits), "]")
        lab3 <- paste0("Q3 (", fmt_cut(b, digits), "–", fmt_cut(c, digits), "]")
        lab4 <- paste0("Q4 > ", fmt_cut(c, digits))
        labs <- c(lab1, lab2, lab3, lab4)
        
        df[[out_var]] <- cut(
            df[[var]],
            breaks = c(-Inf, a, b, c, Inf),
            labels = labs,
            include.lowest = TRUE,
            right = TRUE
        )
        df[[out_var]] <- factor(df[[out_var]], levels = labs, ordered = TRUE)
        df
    }
    
    # ----------------------------- 9) TESTS + TABLE -----------------------------
    test_categorical <- function(var, group) {
        tbl <- table(var, group)
        if (any(dim(tbl) < 2) || sum(rowSums(tbl) > 0) < 2 || sum(colSums(tbl) > 0) < 2) return(NA_real_)
        out <- NA_real_
        tryCatch({
            chi <- suppressWarnings(chisq.test(tbl))
            out <- if (any(chi$expected < 5)) {
                fisher.test(tbl, simulate.p.value = (nrow(tbl) > 2 || ncol(tbl) > 2))$p.value
            } else {
                chi$p.value
            }
        }, error = function(e) out <<- NA_real_)
        out
    }
    
    create_comparison_table <- function(dataset, group_var_name, yes_only_vars = c("MASLD","SLF")) {
        group_var_quo <- rlang::sym(group_var_name)
        
        vars_to_analyze_all <- c(
            "MASLD","SLF",
            "Alb_NLR",
            "NLR","Albumin_g_dL","Albumin_g_L",
            "RIDAGEYR","BMXBMI","Waist_Hip_Ratio",
            "Neutrophil_percent","Lymphocyte_percent",
            "Neutrophil_count","Lymphocyte_count",
            "LBDTCSI","LBDTRSI",
            "LBXSATSI","LBXSASSI",
            "Gender","Race_Ethnicity","Education_Level","Marital_Status",
            "Smoking_Status","Diabetes_Mellitus","Hypertension","Chronic_Kidney_Disease",
            "Cancer","Cardiovascular_Disease"
        )
        vars_to_analyze_all <- intersect(vars_to_analyze_all, names(dataset))
        
        categorical_vars <- vars_to_analyze_all[sapply(dataset[vars_to_analyze_all], function(col) is.factor(col) || is.character(col))]
        numerical_vars <- setdiff(vars_to_analyze_all, categorical_vars)
        
        grp_vec <- pull(dataset, !!group_var_quo)
        if (all(is.na(grp_vec))) return(NULL)
        group_levels <- if (is.factor(grp_vec)) levels(grp_vec) else sort(unique(na.omit(grp_vec)))
        if (length(group_levels) < 2) return(NULL)
        
        group_counts <- table(grp_vec)
        
        header <- c(
            "Variable","join_key","General",
            sapply(group_levels, function(g) paste0(g, " (n=", as.integer(group_counts[as.character(g)]), ")")),
            "p_value"
        )
        results_table <- data.frame(matrix(ncol = length(header), nrow = 0))
        names(results_table) <- header
        
        test_normality <- function(var) {
            v <- na.omit(var)
            if (length(v) >= 3 && length(v) <= 5000) shapiro.test(v)$p.value else NA
        }
        
        always_keep <- c("MASLD","SLF")  # FIX: do not drop these even if >50% missing
        
        for (var_name in c(numerical_vars, categorical_vars)) {
            if (!(var_name %in% c("LBDTRSI", always_keep)) && mean(is.na(dataset[[var_name]])) > 0.5) next
            
            if (var_name %in% numerical_vars) {
                norm_p <- tapply(dataset[[var_name]], grp_vec, test_normality)
                is_normal <- all(!is.na(norm_p) & norm_p > 0.05)
                
                p_value <- tryCatch(
                    if (is_normal) summary(aov(dataset[[var_name]] ~ grp_vec))[[1]][["Pr(>F)"]][1]
                    else kruskal.test(dataset[[var_name]] ~ grp_vec)$p.value,
                    error = function(e) NA_real_
                )
                
                format_stat <- function(x) {
                    if (all(is.na(x))) ""
                    else if (is_normal) paste0(
                        format(round(mean(x, na.rm = TRUE), 2), nsmall = 2), " \u00B1 ",
                        format(round(sd(x, na.rm = TRUE), 2), nsmall = 2)
                    ) else paste0(
                        format(round(median(x, na.rm = TRUE), 2), nsmall = 2), " (",
                        format(round(quantile(x, 0.25, na.rm = TRUE), 2), nsmall = 2), "-",
                        format(round(quantile(x, 0.75, na.rm = TRUE), 2), nsmall = 2), ")"
                    )
                }
                
                new_row <- list(Variable = var_name, join_key = var_name, General = format_stat(dataset[[var_name]]))
                group_stats <- tapply(dataset[[var_name]], grp_vec, format_stat)
                
                for (g in group_levels) {
                    nm <- as.character(g)
                    val <- group_stats[[nm]]
                    if (is.null(val) || is.na(val)) val <- ""
                    new_row[[paste0(nm, " (n=", as.integer(group_counts[nm]), ")")]] <- val
                }
                
                new_row$p_value <- ifelse(is.na(p_value), "NA",
                                          ifelse(p_value < 0.001, "<0.001", format(round(p_value, 3), nsmall = 3)))
                results_table <- rbind(results_table, as.data.frame(new_row, stringsAsFactors = FALSE))
                
            } else {
                p_value <- test_categorical(dataset[[var_name]], grp_vec)
                
                header_row <- as.list(setNames(rep("", ncol(results_table)), names(results_table)))
                header_row$Variable <- paste0(var_name, ", n (%)")
                header_row$join_key <- var_name
                header_row$p_value <- ifelse(is.na(p_value), "NA",
                                             ifelse(p_value < 0.001, "<0.001", format(round(p_value, 3), nsmall = 3)))
                results_table <- rbind(results_table, as.data.frame(header_row, stringsAsFactors = FALSE))
                
                levels_var <- levels(dataset[[var_name]])
                if (is.null(levels_var)) levels_var <- sort(unique(na.omit(dataset[[var_name]])))
                
                # YES-only for selected binaries
                if (var_name %in% yes_only_vars && length(levels_var) >= 2 && any(tolower(levels_var) == "yes")) {
                    levels_to_print <- levels_var[tolower(levels_var) == "yes"]
                } else {
                    levels_to_print <- levels_var
                }
                
                for (lv in levels_to_print) {
                    if (sum(dataset[[var_name]] == lv, na.rm = TRUE) == 0) next
                    
                    level_row <- list(Variable = paste0("  ", lv))
                    level_row$join_key <- paste0(var_name, "___", lv)
                    
                    cnt_g <- sum(!is.na(dataset[[var_name]]) & dataset[[var_name]] == lv)
                    den_g <- sum(!is.na(dataset[[var_name]]))
                    level_row$General <- ifelse(den_g > 0,
                                                paste0(cnt_g, " (", format(round(100 * cnt_g / den_g, 1), nsmall = 1), ")"),
                                                "")
                    
                    for (g in group_levels) {
                        nm <- as.character(g)
                        group_data <- dataset[grp_vec == nm & !is.na(grp_vec), ]
                        
                        den <- sum(!is.na(group_data[[var_name]]))
                        if (den == 0) {
                            level_row[[paste0(nm, " (n=", nrow(group_data), ")")]] <- ""
                        } else {
                            cnt <- sum(!is.na(group_data[[var_name]]) & group_data[[var_name]] == lv)
                            pct <- 100 * cnt / den
                            level_row[[paste0(nm, " (n=", nrow(group_data), ")")]] <- paste0(cnt, " (", format(round(pct, 1), nsmall = 1), ")")
                        }
                    }
                    
                    level_row$p_value <- ""
                    results_table <- rbind(results_table, as.data.frame(level_row, stringsAsFactors = FALSE))
                }
            }
        }
        
        results_table[is.na(results_table)] <- ""
        
        # Pretty labels
        results_table$Variable <- dplyr::case_when(
            results_table$Variable == "MASLD, n (%)" ~ "MASLD, n (%)",
            results_table$Variable == "SLF, n (%)"   ~ "SLF (within MASLD with valid LUXSMED), n (%)",
            TRUE ~ results_table$Variable
        )
        
        results_table
    }
    
    export_table_word <- function(dataset, group_var_name, caption, filename) {
        tab <- create_comparison_table(dataset, group_var_name)
        if (is.null(tab)) {
            cat("Không đủ nhóm để tạo bảng cho", group_var_name, "\n")
            return(invisible(NULL))
        }
        
        ft <- flextable(tab %>% select(-join_key)) %>%
            autofit() %>%
            theme_booktabs() %>%
            bold(i = ~ !grepl("^  ", Variable), j = "Variable", bold = TRUE) %>%
            align(j = 1, align = "left") %>%
            fix_border_issues() %>%
            set_caption(caption)
        
        doc <- read_docx()
        doc <- body_add_flextable(doc, value = ft)
        print(doc, target = filename)
        cat("Đã lưu Word:", normalizePath(filename), "\n")
    }
    
    # ----------------------------- 10) MAIN FLOW -----------------------------
    if (!exists("data")) stop("Thiếu object `data` trong môi trường.")
    
    N0 <- nrow(data)
    
    ex <- apply_exclusion_criteria(data)
    data_cleaned <- ex$data
    
    cat("BẢNG LOG LOẠI TRỪ (theo bước):\n")
    print(ex$log, row.names = FALSE)
    cat("----------------------------------------------------------\n\n")
    
    data_featured <- create_all_features(data_cleaned)
    
    cat("--- BÁO CÁO GIÁ TRỊ THIẾU SAU XỬ LÝ ---\n")
    missing_report <- generate_missing_report(data_featured)
    print(missing_report, row.names = FALSE)
    cat("----------------------------------------------------------\n\n")
    
    # Analytic set:
    # - Alb_NLR available (quartiles)
    # - masld_group defined
    analytic <- data_featured %>%
        mutate(
            MASLD = factor(if_else(masld_group == "MASLD", "Yes", "No"),
                           levels = c("No", "Yes")),
            SLF = factor(case_when(
                masld_group != "MASLD" ~ NA_character_,
                is.na(slf_group) ~ NA_character_,
                slf_group == "SLF" ~ "Yes",
                slf_group == "non-SLF" ~ "No",
                TRUE ~ NA_character_
            ), levels = c("No", "Yes"))
        ) %>%
        filter(!is.na(Alb_NLR), !is.na(masld_group))
    
    cat("TÓM TẮT SỐ LƯỢNG:\n")
    cat("N ban đầu:", N0, "\n")
    cat("N sau loại trừ:", nrow(data_cleaned), "\n")
    cat("N analytic (Alb/NLR + MASLD-defined):", nrow(analytic), "\n")
    cat("-> Loại do Alb/NLR NA:", sum(is.na(data_featured$Alb_NLR)), "\n")
    cat("-> Loại do masld_group NA:", sum(is.na(data_featured$masld_group)), "\n")
    cat("N MASLD trong analytic:", sum(analytic$MASLD == "Yes", na.rm = TRUE), "\n")
    cat("N SLF-defined trong MASLD (analytic):",
        sum(analytic$MASLD == "Yes" & !is.na(analytic$SLF), na.rm = TRUE), "\n\n")
    
    # Quartiles based on analytic set
    cuts <- compute_albnlr_cuts(analytic, var = "Alb_NLR")
    cat("Alb/NLR cutpoints (a, b, c) từ ANALYTIC:", paste(fmt_cut(cuts, 3), collapse = ", "), "\n\n")
    
    analytic <- add_albnlr_quartiles(analytic, cuts = cuts, var = "Alb_NLR", out_var = "AlbNLR_Q", digits = 3)
    
    cat("Phân bố Alb/NLR quartiles:\n")
    print(table(analytic$AlbNLR_Q, useNA = "ifany"))
    cat("\n")
    
    # Export ONE Word table
    export_table_word(
        dataset = analytic,
        group_var_name = "AlbNLR_Q",
        caption = paste0(
            "Bảng: Đặc điểm theo tứ phân vị Alb/NLR (Q1–Q4). ",
            "MASLD tính trên toàn bộ analytic; SLF chỉ tính trong MASLD có LUXSMED hợp lệ (mẫu số thay đổi theo Q)."
        ),
        filename = "Analysis_AlbNLR_Quartiles_OneTable_with_MASLD_SLF.docx"
    )
    
    cat("\n--- TOÀN BỘ QUÁ TRÌNH ĐÃ HOÀN TẤT ---\n")
    
    ```
    
- **Correlation of the Albumin/NLR with MASLD and SLF**
    
    ```r
    # ==============================================================================
    # END-TO-END: Logistic regression (Alb/NLR, Albumin, Neutrophil, Lymphocyte)
    # - REMOVE: CURI, hs-CRP (LBXHSCRP), Uric acid (LBXSUA/LBDSUASI), quartiles of CURI
    # - KEEP required blood columns (exact):
    #     LBXLYPCT  % lymphocyte
    #     LBXNEPCT  % neutrophil
    #     LBDLYMNO  lymphocyte (count)
    #     LBDNENO   neutrophil (count)
    #     LBXSAL    Albumin (g/dL)
    #     LBDSALSI  Albumin (g/L)
    # - NO train/validation split
    # - Outcomes:
    #     MASLD (masld_binary) in whole eligible dataset
    #     SLF (slf_binary) within MASLD only
    # - Models:
    #     Model1 unadjusted
    #     Model2 + sociodemographics
    #     Model3 + cardiometabolic comorbidities
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(officer)
        library(flextable)
        library(tidyverse)
        library(gt)
        library(broom)
        library(stringr)
    })
    
    # ----------------------------- 0) HELPERS -----------------------------
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    ensure_cols <- function(df, cols, default = NA) {
        miss <- setdiff(cols, names(df))
        if (length(miss)) for (m in miss) df[[m]] <- default
        df
    }
    
    # ----------------------------- 1) ALCOHOL -----------------------------
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_cols(data, c("ALQ121", "ALQ130"), default = NA)
        data %>%
            mutate(
                ALQ121 = safe_num(ALQ121),
                ALQ130 = safe_num(ALQ130),
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # ----------------------------- 2) EXCLUSION -----------------------------
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ ---\n")
        initial_rows_total <- nrow(data)
        
        # drop fully-empty columns
        data <- data[, colSums(!is.na(data)) > 0, drop = FALSE]
        
        # ensure key cols exist to avoid errors
        data <- ensure_cols(
            data,
            c("RIAGENDR","RIDAGEYR","LUXSMED","LUXCAPM",
              "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO",
              "LBXSAL","LBDSALSI",
              "ALQ121","ALQ130",
              "RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV"),
            default = NA
        )
        
        # 1) Viral hepatitis
        cat("1. Loại bỏ các trường hợp viêm gan virus...\n")
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR", "LBDHCV")
        hep_vars_exist <- intersect(hep_vars, names(data))
        if (length(hep_vars_exist) > 0) {
            before_n <- nrow(data)
            data <- data %>%
                mutate(across(all_of(hep_vars_exist), safe_num)) %>%
                filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
            cat("-> Đã loại bỏ", before_n - nrow(data), "người do có bằng chứng viêm gan virus.\n\n")
        } else {
            cat("-> CẢNH BÁO: Không có biến viêm gan virus. Bỏ qua.\n\n")
        }
        
        # 2) Require main measurements: FibroScan + Albumin + (Neut&Lymph: counts OR %)
        cat("2. Loại bỏ các hàng thiếu dữ liệu ở các biến đo lường chính...\n")
        before_n <- nrow(data)
        
        data <- data %>% drop_na(LUXSMED, LUXCAPM)
        
        data <- data %>%
            mutate(
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                Albumin_g_dL_tmp = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXLYPCT = safe_num(LBXLYPCT),
                has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO),
                has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT)
            ) %>%
            filter(!is.na(Albumin_g_dL_tmp)) %>%
            filter(has_counts | has_pct)
        
        cat("-> Đã loại bỏ", before_n - nrow(data), "hàng do thiếu FibroScan/Albumin/Neut-Lymph.\n\n")
        
        # 3) Alcohol threshold
        cat("3. Tính toán và loại bỏ các trường hợp tiêu thụ rượu bia vượt ngưỡng...\n")
        data <- calculate_alcohol_consumption(data)
        before_n <- nrow(data)
        data <- data %>%
            mutate(RIAGENDR = safe_num(RIAGENDR)) %>%
            filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        cat("-> Đã loại bỏ", before_n - nrow(data), "hàng do rượu bia vượt ngưỡng.\n\n")
        
        # 4) Age >= 18
        cat("4. Loại bỏ các đối tượng dưới 18 tuổi...\n")
        before_n <- nrow(data)
        data <- data %>%
            mutate(RIDAGEYR = safe_num(RIDAGEYR)) %>%
            filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        cat("-> Đã loại bỏ", before_n - nrow(data), "hàng do tuổi < 18.\n\n")
        
        cat("--- HOÀN TẤT LOẠI TRỪ. Tổng cộng đã loại bỏ:",
            initial_rows_total - nrow(data),
            "hàng. Số hàng còn lại:", nrow(data), "---\n\n")
        data
    }
    
    # ----------------------------- 3) FEATURE ENGINEERING -----------------------------
    create_all_features <- function(data) {
        cat("--- BẮT ĐẦU: TẠO BIẾN (FEATURE ENGINEERING) ---\n")
        
        # Ensure optional covariate columns exist
        cols_to_check <- c(
            "LBXSGL","SMQ040","DIQ050","DIQ070","BPQ040A","KIQ022","KIQ025",
            "MCQ220","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F",
            "BPQ090D","LBXSCR","LBDSCRSI",
            "LBXGLU","LBXGH","DIQ010",
            "BPQ020",
            "DMDEDUC2","DMDMARTZ","SMQ020","RIDRETH1",
            "BMXBMI","BMXWAIST","BMXHIP",
            "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
            "LBDTRSI","LBDHDDSI","LBDTCSI"
        )
        for (col in cols_to_check) if (!col %in% names(data)) data[[col]] <- NA
        
        data <- data %>%
            mutate(across(any_of(cols_to_check), ~ .x))  # no-op, just stable
        
        # --- MASLD / SLF definitions (unchanged logic, removed CURI parts) ---
        data_featured <- data %>%
            mutate(
                # numeric coercions for key vars
                RIAGENDR = safe_num(RIAGENDR),
                RIDAGEYR = safe_num(RIDAGEYR),
                LUXCAPM  = safe_num(LUXCAPM),
                LUXSMED  = safe_num(LUXSMED),
                BMXBMI   = safe_num(BMXBMI),
                BMXWAIST = safe_num(BMXWAIST),
                BMXHIP   = safe_num(BMXHIP),
                LBXGLU   = safe_num(LBXGLU),
                LBXGH    = safe_num(LBXGH),
                BPXOSY1  = safe_num(BPXOSY1), BPXODI1 = safe_num(BPXODI1),
                BPXOSY2  = safe_num(BPXOSY2), BPXODI2 = safe_num(BPXODI2),
                BPXOSY3  = safe_num(BPXOSY3), BPXODI3 = safe_num(BPXODI3),
                BPQ040A  = safe_num(BPQ040A),
                BPQ020   = safe_num(BPQ020),
                BPQ090D  = safe_num(BPQ090D),
                LBDTRSI  = safe_num(LBDTRSI),
                LBDHDDSI = safe_num(LBDHDDSI),
                LBDTCSI  = safe_num(LBDTCSI),
                DIQ010   = safe_num(DIQ010),
                DIQ050   = safe_num(DIQ050),
                DIQ070   = safe_num(DIQ070),
                SMQ020   = safe_num(SMQ020),
                SMQ040   = safe_num(SMQ040),
                DMDEDUC2 = safe_num(DMDEDUC2),
                DMDMARTZ = safe_num(DMDMARTZ),
                RIDRETH1 = safe_num(RIDRETH1),
                KIQ022   = safe_num(KIQ022),
                KIQ025   = safe_num(KIQ025),
                MCQ220   = safe_num(MCQ220),
                MCQ160B  = safe_num(MCQ160B),
                MCQ160C  = safe_num(MCQ160C),
                MCQ160D  = safe_num(MCQ160D),
                MCQ160E  = safe_num(MCQ160E),
                LBXSCR   = safe_num(LBXSCR),
                LBDSCRSI = safe_num(LBDSCRSI)
            ) %>%
            mutate(
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                
                risk_bmi_waist = if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                    NA,
                    (BMXBMI >= 25) |
                        (RIAGENDR == 1 & BMXWAIST >= 94) |
                        (RIAGENDR == 2 & BMXWAIST >= 80)
                ),
                
                risk_glucose_diabetes = if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
                        (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                
                any_bp_measurement_high =
                    (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                
                risk_blood_pressure = if_else(
                    is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
                    NA,
                    any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                
                risk_triglycerides = if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                risk_low_hdl = if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))),
                    NA,
                    ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                )
            ) %>%
            rowwise() %>%
            mutate(
                num_cardiometabolic_risks = sum(
                    c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
                    na.rm = TRUE
                )
            ) %>%
            ungroup() %>%
            mutate(
                has_one_plus_cardiometabolic_risk = if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA,
                    num_cardiometabolic_risks > 0
                ),
                
                is_light_drinker_final = case_when(
                    is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
                    (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                
                masld_group = case_when(
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
                    TRUE ~ NA_character_
                ),
                masld_binary = if_else(masld_group == "MASLD", 1L, 0L, missing = NA_integer_),
                
                slf_group = case_when(
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ "SLF",
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED < 8.0  ~ "non-SLF",
                    TRUE ~ NA_character_
                ),
                slf_binary = if_else(slf_group == "SLF", 1L, 0L, missing = NA_integer_)
            )
        
        # --- Blood indices: Albumin, Neutrophil, Lymphocyte, NLR, Alb/NLR (required columns only) ---
        data_featured <- data_featured %>%
            mutate(
                LBXLYPCT = safe_num(LBXLYPCT),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBDNENO  = safe_num(LBDNENO),
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                
                Albumin_g_dL = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                Albumin_g_L = case_when(
                    !is.na(LBDSALSI) ~ LBDSALSI,
                    is.na(LBDSALSI) & !is.na(LBXSAL) ~ LBXSAL * 10,
                    TRUE ~ NA_real_
                ),
                
                Neutrophil_percent = LBXNEPCT,
                Lymphocyte_percent = LBXLYPCT,
                Neutrophil_count   = LBDNENO,
                Lymphocyte_count   = LBDLYMNO,
                
                NLR = case_when(
                    !is.na(Neutrophil_count) & !is.na(Lymphocyte_count) & Lymphocyte_count > 0 ~
                        Neutrophil_count / Lymphocyte_count,
                    (is.na(Neutrophil_count) | is.na(Lymphocyte_count)) &
                        !is.na(Neutrophil_percent) & !is.na(Lymphocyte_percent) & Lymphocyte_percent > 0 ~
                        Neutrophil_percent / Lymphocyte_percent,
                    TRUE ~ NA_real_
                ),
                Alb_NLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0, Albumin_g_dL / NLR, NA_real_)
            )
        
        # --- Covariates for regression (same structure as your code, renamed) ---
        data_featured <- data_featured %>%
            mutate(
                age = RIDAGEYR,
                gender = factor(RIAGENDR, labels = c("Male", "Female")),
                race = factor(case_when(
                    RIDRETH1 == 1 ~ "Mexican American",
                    RIDRETH1 == 2 ~ "Other Hispanic",
                    RIDRETH1 == 3 ~ "Non-Hispanic White",
                    RIDRETH1 == 4 ~ "Non-Hispanic Black",
                    RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
                    TRUE ~ "Other Race/Multi-Racial"
                )),
                marital_status = factor(case_when(
                    DMDMARTZ %in% c(1, 6) ~ "Married/Living with partner",
                    DMDMARTZ %in% c(2, 3, 4) ~ "Widowed/Divorced/Separated",
                    DMDMARTZ == 5 ~ "Never married",
                    TRUE ~ NA_character_
                )),
                education = factor(case_when(
                    DMDEDUC2 == 1 ~ "Less than 9th grade",
                    DMDEDUC2 == 2 ~ "9-11th grade",
                    DMDEDUC2 == 3 ~ "High school graduate",
                    DMDEDUC2 == 4 ~ "Some college or AA degree",
                    DMDEDUC2 == 5 ~ "College graduate or above",
                    TRUE ~ NA_character_
                )),
                smoking_status = factor(case_when(
                    SMQ020 == 2 ~ "Non-smoker",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
                    TRUE ~ NA_character_
                ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
                
                creatinine_mg_dl = case_when(
                    !is.na(LBXSCR) ~ LBXSCR,
                    is.na(LBXSCR) & !is.na(LBDSCRSI) ~ LBDSCRSI / 88.4,
                    TRUE ~ NA_real_
                ),
                eGFR = if_else(
                    !is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                    175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                        if_else(RIAGENDR == 2, 0.742, 1) * if_else(RIDRETH1 == 4, 1.212, 1),
                    NA_real_
                ),
                
                obesity = factor(if_else(BMXBMI >= 30, "Yes", "No", missing = "No"), levels = c("No","Yes")),
                diabetes_mellitus = factor(
                    if_else(DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                                (!is.na(LBXGH) & LBXGH >= 6.5) |
                                (!is.na(LBXSGL) & LBXSGL >= 126),
                            "Yes","No", missing = "No"),
                    levels = c("No","Yes")
                ),
                hypertension = factor(if_else(BPQ020 == 1 | BPQ040A == 1, "Yes","No", missing = "No"), levels = c("No","Yes")),
                chronic_kidney_disease = factor(if_else((KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60)), "Yes","No", missing = "No"), levels = c("No","Yes")),
                cancer = factor(if_else(MCQ220 == 1, "Yes","No", missing = "No"), levels = c("No","Yes")),
                cardiovascular_disease = factor(if_else(MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1, "Yes","No", missing = "No"), levels = c("No","Yes"))
            )
        
        cat("--- HOÀN TẤT TẠO BIẾN. ---\n\n")
        data_featured
    }
    
    # ----------------------------- 4) REGRESSION (Alb/NLR etc.) -----------------------------
    run_regression_analysis <- function(data) {
        
        analyze_outcome <- function(outcome_var, data, predictor_spec) {
            
            covariate_sets <- list(
                Model1 = c(),
                Model2 = c("age", "gender", "race", "marital_status", "education"),
                Model3 = c("age", "gender", "race", "marital_status", "education",
                           "smoking_status", "obesity", "diabetes_mellitus", "hypertension",
                           "cancer", "chronic_kidney_disease", "cardiovascular_disease")
            )
            
            format_results <- function(model, term_name) {
                tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
                    filter(str_detect(term, term_name)) %>%
                    mutate(
                        OR_95CI = paste0(sprintf("%.2f", estimate), " (", sprintf("%.2f", conf.low), ", ", sprintf("%.2f", conf.high), ")"),
                        p.value_formatted = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
                    ) %>%
                    dplyr::select(term, OR_95CI, p.value_formatted)
            }
            
            all_results <- list()
            
            for (model_name in names(covariate_sets)) {
                
                predictor <- predictor_spec$predictor
                predictor_label <- predictor_spec$label
                predictor_type <- predictor_spec$type  # "continuous" or "quartile"
                
                res_main <- tibble()
                res_quart <- tibble()
                res_trend <- tibble()
                
                if (predictor_type == "continuous") {
                    tryCatch({
                        covs <- covariate_sets[[model_name]]
                        data_model <- data %>% select(any_of(c(outcome_var, predictor, covs))) %>% na.omit()
                        
                        if (nrow(data_model) < 20) stop("Too few complete cases.")
                        if (n_distinct(data_model[[outcome_var]]) < 2) stop("Outcome variable has fewer than 2 levels.")
                        
                        valid_covs <- covs[sapply(covs, function(cov)
                            if (!cov %in% names(data_model)) FALSE
                            else !is.factor(data_model[[cov]]) || n_distinct(data_model[[cov]]) >= 2
                        )]
                        formula_part <- if (length(valid_covs) > 0) paste(valid_covs, collapse = " + ") else "1"
                        fml <- as.formula(paste(outcome_var, "~", predictor, "+", formula_part))
                        
                        model <- glm(fml, data = data_model, family = binomial())
                        res_main <- format_results(model, predictor)
                    }, error = function(e) {
                        cat("--> ERROR (", predictor_label, " ", model_name, "): ", e$message, "\n", sep="")
                    })
                    
                    all_results[[model_name]] <- bind_rows(
                        tibble(Category = predictor_label,
                               OR = if (nrow(res_main) > 0) res_main$OR_95CI[1] else "Lỗi",
                               P  = if (nrow(res_main) > 0) res_main$p.value_formatted[1] else "")
                    )
                    
                } else if (predictor_type == "quartile") {
                    # quartile with Q1 ref + trend
                    predictor_quart <- predictor_spec$predictor
                    predictor_num   <- predictor_spec$trend_numeric  # numeric-coded quartile variable name
                    tryCatch({
                        covs <- covariate_sets[[model_name]]
                        data_quart <- data %>%
                            select(any_of(c(outcome_var, predictor_quart, predictor_num, covs))) %>%
                            na.omit()
                        
                        if (nrow(data_quart) < 40) stop("Too few complete cases.")
                        if (n_distinct(data_quart[[outcome_var]]) < 2) stop("Outcome variable has fewer than 2 levels.")
                        if (n_distinct(data_quart[[predictor_quart]]) < 2) stop("Quartile variable has fewer than 2 levels.")
                        if (!"Q1" %in% levels(data_quart[[predictor_quart]])) stop("Reference level 'Q1' not found.")
                        
                        valid_covs <- covs[sapply(covs, function(cov)
                            if (!cov %in% names(data_quart)) FALSE
                            else !is.factor(data_quart[[cov]]) || n_distinct(data_quart[[cov]]) >= 2
                        )]
                        formula_part <- if (length(valid_covs) > 0) paste(valid_covs, collapse = " + ") else "1"
                        
                        data_quart[[predictor_quart]] <- relevel(data_quart[[predictor_quart]], ref = "Q1")
                        
                        fml_quart <- as.formula(paste(outcome_var, "~", predictor_quart, "+", formula_part))
                        m_quart <- glm(fml_quart, data = data_quart, family = binomial())
                        res_quart <- format_results(m_quart, predictor_quart)
                        
                        fml_trend <- as.formula(paste(outcome_var, "~", predictor_num, "+", formula_part))
                        m_trend <- glm(fml_trend, data = data_quart, family = binomial())
                        res_trend <- format_results(m_trend, predictor_num)
                        
                    }, error = function(e) {
                        cat("--> ERROR (", predictor_label, " quartile ", model_name, "): ", e$message, "\n", sep="")
                    })
                    
                    # res_quart typically has 3 rows for Q2/Q3/Q4 vs Q1
                    all_results[[model_name]] <- bind_rows(
                        tibble(Category = paste0(predictor_label, " quartiles"), OR = "", P = ""),
                        tibble(Category = "Q1", OR = "Reference", P = ""),
                        tibble(Category = "Q2",
                               OR = if (nrow(res_quart) >= 1) res_quart$OR_95CI[1] else NA_character_,
                               P  = if (nrow(res_quart) >= 1) res_quart$p.value_formatted[1] else NA_character_),
                        tibble(Category = "Q3",
                               OR = if (nrow(res_quart) >= 2) res_quart$OR_95CI[2] else NA_character_,
                               P  = if (nrow(res_quart) >= 2) res_quart$p.value_formatted[2] else NA_character_),
                        tibble(Category = "Q4",
                               OR = if (nrow(res_quart) >= 3) res_quart$OR_95CI[3] else NA_character_,
                               P  = if (nrow(res_quart) >= 3) res_quart$p.value_formatted[3] else NA_character_),
                        tibble(Category = "P for trend", OR = "", P = if (nrow(res_trend) >= 1) res_trend$p.value_formatted[1] else NA_character_)
                    )
                }
            }
            
            # merge columns per model
            final_df <- Reduce(function(x, y) left_join(x, y, by = "Category"),
                               lapply(names(all_results), function(nm) {
                                   df <- all_results[[nm]]
                                   colnames(df) <- c("Category", paste0(nm, "_OR"), paste0(nm, "_P"))
                                   df
                               }))
            final_df
        }
        
        # Predictor blocks:
        # 1) Alb_NLR continuous (primary)
        # 2) Albumin_g_dL continuous
        # 3) NLR continuous
        # 4) NLR quartiles with trend
        # 5) Alb_NLR quartiles with trend
        make_quartiles <- function(x, labels = c("Q1","Q2","Q3","Q4")) {
            qs <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7)
            if (any(!is.finite(qs)) || length(unique(qs)) < 5) return(rep(NA_character_, length(x)))
            as.character(cut(x, breaks = qs, include.lowest = TRUE, labels = labels))
        }
        
        data <- data %>%
            mutate(
                Alb_NLR_q = factor(make_quartiles(Alb_NLR), levels = c("Q1","Q2","Q3","Q4")),
                NLR_q     = factor(make_quartiles(NLR),     levels = c("Q1","Q2","Q3","Q4")),
                Alb_NLR_q_num = if_else(!is.na(Alb_NLR_q), as.numeric(Alb_NLR_q), NA_real_),
                NLR_q_num     = if_else(!is.na(NLR_q),     as.numeric(NLR_q),     NA_real_)
            )
        
        # --- MASLD models (whole eligible dataset) ---
        cat("--- Bắt đầu phân tích hồi quy cho MASLD ---\n")
        masld_header <- tibble(Category = "MASLD", Model1_OR = "", Model1_P = "", Model2_OR = "", Model2_P = "", Model3_OR = "", Model3_P = "")
        
        # continuous predictors
        masld_AlbNLR <- analyze_outcome("masld_binary", data, list(predictor="Alb_NLR", label="Alb/NLR (continuous)", type="continuous"))
        masld_Alb    <- analyze_outcome("masld_binary", data, list(predictor="Albumin_g_dL", label="Albumin (g/dL, continuous)", type="continuous"))
        masld_NLR    <- analyze_outcome("masld_binary", data, list(predictor="NLR", label="NLR (continuous)", type="continuous"))
        
        # quartiles
        masld_NLR_q  <- analyze_outcome("masld_binary", data, list(predictor="NLR_q", label="NLR", type="quartile", trend_numeric="NLR_q_num"))
        masld_AlbNLR_q <- analyze_outcome("masld_binary", data, list(predictor="Alb_NLR_q", label="Alb/NLR", type="quartile", trend_numeric="Alb_NLR_q_num"))
        
        masld_block <- bind_rows(masld_header, masld_AlbNLR, masld_Alb, masld_NLR, masld_NLR_q, masld_AlbNLR_q)
        
        # --- SLF models (MASLD only) ---
        cat("--- Bắt đầu phân tích hồi quy cho SLF (chỉ trong nhóm MASLD) ---\n")
        data_slf_only <- data %>% filter(masld_group == "MASLD")
        
        slf_header <- tibble(Category = "Significant clinical fibrosis", Model1_OR = "", Model1_P = "", Model2_OR = "", Model2_P = "", Model3_OR = "", Model3_P = "")
        
        slf_AlbNLR <- analyze_outcome("slf_binary", data_slf_only, list(predictor="Alb_NLR", label="Alb/NLR (continuous)", type="continuous"))
        slf_Alb    <- analyze_outcome("slf_binary", data_slf_only, list(predictor="Albumin_g_dL", label="Albumin (g/dL, continuous)", type="continuous"))
        slf_NLR    <- analyze_outcome("slf_binary", data_slf_only, list(predictor="NLR", label="NLR (continuous)", type="continuous"))
        slf_NLR_q  <- analyze_outcome("slf_binary", data_slf_only, list(predictor="NLR_q", label="NLR", type="quartile", trend_numeric="NLR_q_num"))
        slf_AlbNLR_q <- analyze_outcome("slf_binary", data_slf_only, list(predictor="Alb_NLR_q", label="Alb/NLR", type="quartile", trend_numeric="Alb_NLR_q_num"))
        
        slf_block <- bind_rows(slf_header, slf_AlbNLR, slf_Alb, slf_NLR, slf_NLR_q, slf_AlbNLR_q)
        
        # combine and format to GT table
        final_table_data <- bind_rows(masld_block, slf_block) %>%
            # standardize colnames expected by gt
            mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))
        
        # Ensure exact column order (Category + per model)
        # analyze_outcome returns Category + Model1_OR/Model1_P ...
        # If any missing columns due to merge edge cases, create them.
        needed <- c("Category","Model1_OR","Model1_P","Model2_OR","Model2_P","Model3_OR","Model3_P")
        for (nm in setdiff(needed, names(final_table_data))) final_table_data[[nm]] <- ""
        
        final_table_data <- final_table_data %>% select(all_of(needed))
        colnames(final_table_data) <- c("Category","OR_1","P_1","OR_2","P_2","OR_3","P_3")
        
        final_gt_table <- final_table_data %>%
            gt() %>%
            tab_header(title = "Logistic Regression of Alb/NLR, Albumin, Neutrophil, Lymphocyte-derived Indices on Outcomes") %>%
            tab_spanner(label = "Model 1", columns = c(OR_1, P_1)) %>%
            tab_spanner(label = "Model 2", columns = c(OR_2, P_2)) %>%
            tab_spanner(label = "Model 3", columns = c(OR_3, P_3)) %>%
            cols_label(
                Category = "",
                OR_1 = "OR (95% CI)", P_1 = "P-value",
                OR_2 = "OR (95% CI)", P_2 = "P-value",
                OR_3 = "OR (95% CI)", P_3 = "P-value"
            ) %>%
            tab_style(
                style = cell_text(weight = "bold"),
                locations = cells_body(columns = Category, rows = Category %in% c("MASLD", "Significant clinical fibrosis"))
            ) %>%
            cols_align(align = "left", columns = Category) %>%
            cols_align(align = "center", columns = -Category) %>%
            tab_source_note(md("**Model 1**: Unadjusted.")) %>%
            tab_source_note(md("**Model 2**: Adjusted for age, gender, race, marital status, and education.")) %>%
            tab_source_note(md("**Model 3**: Adjusted for Model 2 variables + smoking status, obesity, diabetes mellitus, hypertension, cancer, chronic kidney disease, and cardiovascular disease."))
        
        final_gt_table
    }
    
    # ----------------------------- 5) MAIN FLOW (NO TRAIN/VAL) -----------------------------
    if (!exists("data")) stop("Thiếu object `data` trong môi trường. Gán: data <- MASLD")
    
    data_cleaned  <- apply_exclusion_criteria(data)
    data_featured <- create_all_features(data_cleaned)
    
    # Keep rows with definable outcomes + key predictors (Alb_NLR / Albumin / NLR)
    data_complete <- data_featured %>%
        filter(!is.na(masld_group)) %>%
        # need at least one of the key predictors present; Alb_NLR preferred
        filter(!is.na(Alb_NLR) | (!is.na(Albumin_g_dL) & !is.na(NLR)))
    
    cat("Tổng số quan sát hợp lệ cho hồi quy:", nrow(data_complete), "\n")
    
    final_table <- run_regression_analysis(data_complete)
    
    cat("\n--- TOÀN BỘ QUÁ TRÌNH PHÂN TÍCH ĐÃ HOÀN TẤT! ---\n")
    print(final_table)
    
    ```
    
- **GAM**
    
    ```r
    # ==============================================================================
    # END-TO-END: GAM (no train/validation) + AUC
    # - Apply exclusion criteria
    # - Create MASLD/SLF outcomes
    # - REMOVE: CURI, hs-CRP, uric acid everywhere
    # - Use ONLY specified blood columns for ANLR construction
    # - Predictor:
    #      ANLR = Albumin_g_dL / NLR
    # - Fit GAM: outcome ~ s(ANLR, bs="cs") on full eligible dataset
    # - Plot with CI + rugs (controls bottom, cases top)
    # - Compute AUC from predicted probabilities (pROC-safe)
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(mgcv)
        library(ggplot2)
        library(pROC)
    })
    
    # ----------------------------
    # Helpers
    # ----------------------------
    ensure_vars <- function(df, vars) {
        miss <- setdiff(vars, names(df))
        if (length(miss) > 0) for (v in miss) df[[v]] <- NA
        df
    }
    
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    # =========================
    # 1) Alcohol consumption
    # =========================
    calculate_alcohol_consumption <- function(data) {
        need <- c("ALQ121", "ALQ130")
        data <- ensure_vars(data, need)
        
        data %>%
            mutate(
                ALQ121 = safe_num(ALQ121),
                ALQ130 = safe_num(ALQ130),
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # =========================
    # 2) Exclusion criteria
    # =========================
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ ---\n")
        initial_rows_total <- nrow(data)
        
        data <- data %>% dplyr::select(where(~ any(!is.na(.x))))
        data <- ensure_vars(data, c("RIAGENDR", "RIDAGEYR", "LUXSMED", "LUXCAPM", "ALQ121", "ALQ130",
                                    "LBXLYPCT", "LBXNEPCT", "LBDLYMNO", "LBDNENO", "LBXSAL", "LBDSALSI"))
        
        # 1) Viral hepatitis
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR", "LBDHCV")
        hep_vars_exist <- intersect(hep_vars, names(data))
        if (length(hep_vars_exist) > 0) {
            data <- data %>%
                mutate(across(all_of(hep_vars_exist), safe_num)) %>%
                filter(rowSums(across(all_of(hep_vars_exist), ~ (.x == 1)), na.rm = TRUE) == 0)
        }
        
        # 2) Drop NA on key measurements
        data <- data %>%
            mutate(
                RIAGENDR = safe_num(RIAGENDR), RIDAGEYR = safe_num(RIDAGEYR),
                LUXSMED  = safe_num(LUXSMED), LUXCAPM  = safe_num(LUXCAPM),
                LBXSAL   = safe_num(LBXSAL), LBDSALSI = safe_num(LBDSALSI),
                LBDNENO  = safe_num(LBDNENO), LBDLYMNO = safe_num(LBDLYMNO),
                LBXNEPCT = safe_num(LBXNEPCT), LBXLYPCT = safe_num(LBXLYPCT),
                Albumin_g_dL_tmp = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO),
                has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT)
            ) %>%
            drop_na(LUXSMED, LUXCAPM) %>%
            filter(!is.na(Albumin_g_dL_tmp)) %>%
            filter(has_counts | has_pct)
        
        # 3) Alcohol
        data <- calculate_alcohol_consumption(data)
        data <- data %>%
            filter(is.na(weekly_alcohol_grams) |
                       (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                       (RIAGENDR == 2 & weekly_alcohol_grams <= 140))
        
        # 4) Age >= 18
        data <- data %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        cat("--- HOÀN TẤT LOẠI TRỪ. Số hàng còn lại:", nrow(data), "---\n\n")
        data
    }
    
    # =========================
    # 3) Feature engineering
    # =========================
    create_features <- function(data) {
        cat("--- BẮT ĐẦU: TẠO BIẾN CHẨN ĐOÁN + ANLR ---\n")
        need <- c("RIAGENDR","RIDAGEYR","BMXBMI","BMXWAIST","LBXGLU","LBXGH","DIQ010","DIQ050","DIQ070",
                  "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3","BPQ040A","BPQ020",
                  "LBDTRSI","BPQ090D","LBDHDDSI","LUXCAPM","LUXSMED","weekly_alcohol_grams",
                  "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI")
        data <- ensure_vars(data, need)
        
        data_processed <- data %>%
            mutate(across(all_of(need), safe_num)) %>%
            mutate(
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                risk_bmi_waist = if_else(is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)), NA,
                                         (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)),
                risk_glucose_diabetes = if_else(is.na(LBXGLU) & is.na(LBXGH) & (is.na(DIQ010) | DIQ010 %in% c(7, 9)), NA,
                                                (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) | (if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
                                                    (if_else(is.na(DIQ010), FALSE, DIQ010 == 1))),
                any_bp_measurement_high = (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                risk_blood_pressure = if_else(is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)), NA,
                                              any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))),
                risk_triglycerides = if_else(is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)), NA,
                                             (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) | (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))),
                risk_low_hdl = if_else(is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))), NA,
                                       ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                                            (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) | (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))),
                num_cardiometabolic_risks = rowSums(cbind(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl), na.rm = TRUE),
                has_one_plus_cardiometabolic_risk = if_else(is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes), NA, num_cardiometabolic_risks > 0),
                is_light_drinker_final = case_when(is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
                                                   (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
                                                   TRUE ~ FALSE),
                masld_group = case_when(is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                                        has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                                        has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == TRUE ~ "MASLD",
                                        TRUE ~ "non-MASLD"),
                masld_binary = if_else(masld_group == "MASLD", 1L, 0L),
                slf_group = case_when(masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ "SLF",
                                      masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED <  8.0 ~ "non-SLF",
                                      TRUE ~ NA_character_),
                slf_binary = if_else(slf_group == "SLF", 1L, 0L),
                
                # === ANLR Calculation ===
                Albumin_g_dL = case_when(!is.na(LBXSAL) ~ LBXSAL, is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10, TRUE ~ NA_real_),
                NLR = case_when(!is.na(LBDNENO) & !is.na(LBDLYMNO) & LBDLYMNO > 0 ~ LBDNENO / LBDLYMNO,
                                (is.na(LBDNENO) | is.na(LBDLYMNO)) & !is.na(LBXNEPCT) & !is.na(LBXLYPCT) & LBXLYPCT > 0 ~ LBXNEPCT / LBXLYPCT,
                                TRUE ~ NA_real_),
                ANLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0, Albumin_g_dL / NLR, NA_real_)
            )
        data_processed
    }
    
    # =========================
    # 4) GAM + plot + AUC
    # =========================
    run_gam_analysis_for_outcome <- function(data, outcome_var, outcome_name, y_limit, y_breaks) {
        cat(paste0("\n--- Phân tích GAM cho ", outcome_name, " ---\n"))
        model_data <- data %>% filter(!is.na(ANLR) & is.finite(ANLR) & !is.na(.data[[outcome_var]]))
        if (nrow(model_data) > 20) {
            cut99 <- quantile(model_data$ANLR, 0.99, na.rm = TRUE, names = FALSE)
            model_data <- model_data %>% filter(ANLR < cut99)
        }
        
        gam_model <- mgcv::gam(as.formula(paste0(outcome_var, " ~ s(ANLR, bs = \"cs\")")), 
                               data = model_data, family = binomial, method = "REML")
        print(summary(gam_model))
        
        # AUC
        probs <- predict(gam_model, type = "response")
        roc_obj <- pROC::roc(response = model_data[[outcome_var]], predictor = probs, quiet = TRUE)
        cat("AUC for ", outcome_name, ": ", round(as.numeric(pROC::auc(roc_obj)), 4), "\n", sep = "")
        
        # Plotting
        plot_data <- data.frame(ANLR = seq(min(model_data$ANLR), max(model_data$ANLR), length.out = 200))
        preds <- predict(gam_model, newdata = plot_data, type = "link", se.fit = TRUE)
        plot_data <- plot_data %>% mutate(prob_fit = plogis(preds$fit),
                                          prob_lower = plogis(preds$fit - 1.96 * preds$se.fit),
                                          prob_upper = plogis(preds$fit + 1.96 * preds$se.fit))
        
        ggplot(plot_data, aes(x = ANLR)) +
            geom_line(aes(y = prob_lower), color = "gray50", linetype = "dotted") +
            geom_line(aes(y = prob_upper), color = "gray50", linetype = "dotted") +
            geom_line(aes(y = prob_fit), color = "red", linewidth = 1.2) +
            geom_rug(data = filter(model_data, .data[[outcome_var]] == 0), aes(x = ANLR, color = "Control"), sides = "b", alpha = 0.3) +
            geom_rug(data = filter(model_data, .data[[outcome_var]] == 1), aes(x = ANLR, color = "Case"), sides = "t", alpha = 0.3) +
            scale_y_continuous(limits = y_limit, breaks = y_breaks) +
            scale_color_manual(values = c("Control" = "grey80", "Case" = "#003366")) +
            labs(x = "ANLR", y = paste("Probability of", outcome_name), color = "Status") +
            theme_classic(base_size = 14)
    }
    
    # ============================================================
    # MAIN
    # ============================================================
    stopifnot(exists("data"))
    data_cleaned  <- apply_exclusion_criteria(data)
    data_featured <- create_features(data_cleaned)
    
    gam_plot_masld <- run_gam_analysis_for_outcome(
        data = data_featured, outcome_var = "masld_binary", outcome_name = "MASLD",
        y_limit = c(0, 1.0), y_breaks = seq(0, 1, 0.2)
    )
    
    gam_plot_slf <- run_gam_analysis_for_outcome(
        data = data_featured %>% filter(masld_group == "MASLD"), 
        outcome_var = "slf_binary", outcome_name = "SLF",
        y_limit = c(0, 0.5), y_breaks = seq(0, 0.5, 0.1)
    )
    
    print(gam_plot_masld)
    print(gam_plot_slf)
    ```
    
- **Youden**
    
    ```
    --- Phân tích GAM và tìm điểm cắt cho: A: MASLD (Full cohort) ---
    AUC (GAM fitted) = 0.5346
    ==> Breakpoint (Youden) = 0.606
    
    --- Phân tích GAM và tìm điểm cắt cho: B: SLF (MASLD population) ---
    AUC (GAM fitted) = 0.5737
    ==> Breakpoint (Youden) = 5.642
    ```
    
    ```r
    # ==============================================================================
    # END-TO-END: GAM + Youden breakpoint (NO train/validation)
    # - Exclusion criteria uses ONLY Alb/NLR inputs as "key measurements" (NO hs-CRP, NO uric acid, NO WBC)
    # - Replace ALL CURI/hs-CRP/uric acid references with Alb_by_NLR
    # - Use exactly these columns (as you specified):
    #     LBXLYPCT  % lymphocyte
    #     LBXNEPCT  % neutrophil
    #     LBDLYMNO  lymphocyte (absolute)
    #     LBDNENO   neutrophil (absolute)
    #     LBXSAL    Albumin (g/dL)
    #     LBDSALSI  Albumin (g/L)
    # - NLR priority: counts (LBDNENO/LBDLYMNO). Fallback: % (LBXNEPCT/LBXLYPCT) if counts missing
    # - Albumin: priority LBXSAL, fallback LBDSALSI/10
    # - Outcomes:
    #     MASLD: CAP>=263 + >=1 cardiometabolic risk (+ light drinker logic preserved)
    #     SLF: within MASLD only, LSM>=8.0
    # - Plots: 2 panels (A MASLD, B SLF), no split
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(mgcv)
        library(ggplot2)
        library(patchwork)
        library(pROC)
    })
    
    # ----------------------------
    # Helpers
    # ----------------------------
    ensure_vars <- function(df, vars) {
        miss <- setdiff(vars, names(df))
        if (length(miss) > 0) for (v in miss) df[[v]] <- NA
        df
    }
    
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    # =========================
    # 1) Alcohol consumption
    # =========================
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_vars(data, c("ALQ121", "ALQ130"))
        
        data %>%
            mutate(
                ALQ121 = safe_num(ALQ121),
                ALQ130 = safe_num(ALQ130),
                
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # =========================
    # 2) Exclusion criteria (key measurements = Alb/NLR inputs + FibroScan)
    # =========================
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ ---\n")
        initial_rows_total <- nrow(data)
        
        # drop fully-empty columns
        data <- data %>% dplyr::select(where(~ any(!is.na(.x))))
        
        # ensure required columns exist (NO hs-CRP, NO uric acid)
        req <- c(
            "RIAGENDR","RIDAGEYR","LUXSMED","LUXCAPM","ALQ121","ALQ130",
            "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"
        )
        data <- ensure_vars(data, req)
        
        # 1) Viral hepatitis
        cat("1. Loại bỏ các trường hợp viêm gan virus...\n")
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR", "LBDHCV")
        hep_vars_exist <- intersect(hep_vars, names(data))
        
        rows_before <- nrow(data)
        if (length(hep_vars_exist) > 0) {
            data <- data %>%
                mutate(across(all_of(hep_vars_exist), safe_num)) %>%
                filter(rowSums(across(all_of(hep_vars_exist), ~ (.x == 1)), na.rm = TRUE) == 0)
        }
        cat("-> Đã loại bỏ", rows_before - nrow(data), "người do có bằng chứng viêm gan virus.\n\n")
        
        # 2) Key measurements: FibroScan + Albumin + Neut/Lymph (counts OR %)
        cat("2. Loại bỏ các hàng thiếu dữ liệu ở các biến đo lường chính (FibroScan + Alb/NLR inputs)...\n")
        rows_before <- nrow(data)
        
        data <- data %>%
            mutate(
                RIAGENDR = safe_num(RIAGENDR),
                RIDAGEYR = safe_num(RIDAGEYR),
                LUXSMED  = safe_num(LUXSMED),
                LUXCAPM  = safe_num(LUXCAPM),
                
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXLYPCT = safe_num(LBXLYPCT),
                
                Albumin_g_dL_tmp = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO) & LBDLYMNO > 0,
                has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT) & LBXLYPCT > 0
            ) %>%
            drop_na(LUXSMED, LUXCAPM) %>%
            filter(!is.na(Albumin_g_dL_tmp)) %>%
            filter(has_counts | has_pct)
        
        cat("-> Đã loại bỏ", rows_before - nrow(data), "hàng do thiếu FibroScan/Albumin/Neut-Lymph.\n\n")
        
        # 3) Alcohol threshold
        cat("3. Tính toán và loại bỏ các trường hợp tiêu thụ rượu bia vượt ngưỡng...\n")
        data <- calculate_alcohol_consumption(data)
        
        rows_before <- nrow(data)
        data <- data %>%
            filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        cat("-> Đã loại bỏ", rows_before - nrow(data), "hàng do tiêu thụ rượu bia vượt ngưỡng.\n\n")
        
        # 4) Age >= 18
        cat("4. Loại bỏ các đối tượng dưới 18 tuổi...\n")
        rows_before <- nrow(data)
        data <- data %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        cat("-> Đã loại bỏ", rows_before - nrow(data), "hàng do tuổi < 18.\n\n")
        
        cat("--- HOÀN TẤT LOẠI TRỪ. Tổng cộng đã loại bỏ:",
            initial_rows_total - nrow(data),
            "hàng. Số hàng còn lại:", nrow(data), "---\n\n")
        
        data
    }
    
    # =========================
    # 3) Feature engineering (MASLD/SLF + Alb/NLR)
    # =========================
    create_features <- function(data) {
        cat("--- BẮT ĐẦU: TẠO BIẾN CHẨN ĐOÁN (MASLD/SLF) + Alb/NLR ---\n")
        
        need <- c(
            "RIAGENDR","RIDAGEYR",
            "BMXBMI","BMXWAIST","LBXGLU","LBXGH","DIQ010","DIQ050","DIQ070",
            "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3","BPQ040A","BPQ020",
            "LBDTRSI","BPQ090D","LBDHDDSI",
            "LUXCAPM","LUXSMED","weekly_alcohol_grams",
            # exact Alb/NLR inputs
            "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"
        )
        data <- ensure_vars(data, need)
        
        d <- data %>%
            mutate(across(all_of(need), safe_num)) %>%
            mutate(
                # --- MASLD definition ---
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                
                risk_bmi_waist = if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                    NA,
                    (BMXBMI >= 25) |
                        (RIAGENDR == 1 & BMXWAIST >= 94) |
                        (RIAGENDR == 2 & BMXWAIST >= 80)
                ),
                
                risk_glucose_diabetes = if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (if_else(is.na(LBXGH),  FALSE, LBXGH  >= 5.7)) |
                        (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                
                any_bp_measurement_high =
                    (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                
                risk_blood_pressure = if_else(
                    is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
                    NA,
                    any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                
                risk_triglycerides = if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                risk_low_hdl = if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))),
                    NA,
                    ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                num_cardiometabolic_risks = rowSums(
                    cbind(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
                    na.rm = TRUE
                ),
                
                has_one_plus_cardiometabolic_risk = if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA, num_cardiometabolic_risks > 0
                ),
                
                is_light_drinker_final = case_when(
                    is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
                    (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                
                masld_group = case_when(
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
                    TRUE ~ NA_character_
                ),
                masld_binary = if_else(masld_group == "MASLD", 1L, 0L),
                
                # --- SLF only within MASLD ---
                slf_group = case_when(
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ "SLF",
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED <  8.0 ~ "non-SLF",
                    TRUE ~ NA_character_
                ),
                slf_binary = if_else(slf_group == "SLF", 1L, 0L),
                
                # --- Albumin in g/dL and g/L (both kept) ---
                Albumin_g_dL = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                Albumin_g_L = case_when(
                    !is.na(LBDSALSI) ~ LBDSALSI,
                    is.na(LBDSALSI) & !is.na(LBXSAL) ~ LBXSAL * 10,
                    TRUE ~ NA_real_
                ),
                
                # --- Neut/Lymph: counts + % ---
                Neutrophil_abs = LBDNENO,
                Lymphocyte_abs = LBDLYMNO,
                Neutrophil_pct = LBXNEPCT,
                Lymphocyte_pct = LBXLYPCT,
                
                # --- NLR (counts primary; % fallback) ---
                NLR = case_when(
                    !is.na(Neutrophil_abs) & !is.na(Lymphocyte_abs) & Lymphocyte_abs > 0 ~ Neutrophil_abs / Lymphocyte_abs,
                    (is.na(Neutrophil_abs) | is.na(Lymphocyte_abs)) &
                        !is.na(Neutrophil_pct) & !is.na(Lymphocyte_pct) & Lymphocyte_pct > 0 ~ Neutrophil_pct / Lymphocyte_pct,
                    TRUE ~ NA_real_
                ),
                
                Alb_by_NLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                                     Albumin_g_dL / NLR, NA_real_)
            )
        
        cat("--- HOÀN TẤT TẠO BIẾN. ---\n\n")
        d
    }
    
    # =========================
    # 4) Youden cutoff (ROC)
    # =========================
    find_youden_cut <- function(y, x) {
        ok <- !is.na(y) & !is.na(x) & is.finite(x)
        y <- as.integer(y[ok])
        x <- x[ok]
        if (length(y) < 20) return(NA_real_)
        if (length(unique(y)) < 2) return(NA_real_)
        if (!all(unique(y) %in% c(0L, 1L))) return(NA_real_)
        
        roc_obj <- tryCatch(
            pROC::roc(response = y, predictor = x, levels = c(0, 1), direction = "<", quiet = TRUE),
            error = function(e) NULL
        )
        if (is.null(roc_obj)) return(NA_real_)
        
        thr <- tryCatch(
            pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold")[1],
            error = function(e) NA_real_
        )
        as.numeric(thr)
    }
    
    # =========================
    # 5) AUC helper
    # =========================
    safe_auc_pROC <- function(y, p) {
        ok <- !is.na(y) & !is.na(p) & is.finite(p)
        y <- as.integer(y[ok])
        p <- p[ok]
        if (length(y) < 10) return(NA_real_)
        if (length(unique(y)) < 2) return(NA_real_)
        roc_obj <- pROC::roc(response = y, predictor = p, levels = c(0, 1), direction = "<", quiet = TRUE)
        as.numeric(pROC::auc(roc_obj))
    }
    
    # =========================
    # 6) GAM + plot + breakpoint annotation
    # =========================
    run_gam_panel <- function(data, outcome_var, panel_title, y_label) {
        cat("\n---", panel_title, "---\n")
        
        model_data <- data %>%
            filter(!is.na(.data[[outcome_var]]), .data[[outcome_var]] %in% c(0, 1),
                   !is.na(Alb_by_NLR), is.finite(Alb_by_NLR)) %>%
            mutate(outcome = as.integer(.data[[outcome_var]])) %>%
            select(outcome, Alb_by_NLR)
        
        # trim 99th
        if (nrow(model_data) >= 20) {
            cut99 <- quantile(model_data$Alb_by_NLR, 0.99, na.rm = TRUE, names = FALSE)
            model_data <- model_data %>% filter(Alb_by_NLR < cut99)
        }
        
        if (nrow(model_data) < 20 || n_distinct(model_data$outcome) < 2) {
            cat("Không đủ dữ liệu.\n")
            return(ggplot() + theme_void() + labs(title = panel_title))
        }
        
        gam_model <- mgcv::gam(outcome ~ s(Alb_by_NLR, bs = "cs"),
                               data = model_data, family = binomial, method = "REML")
        print(summary(gam_model))
        
        fitted_probs <- predict(gam_model, type = "response")
        auc_val <- safe_auc_pROC(model_data$outcome, fitted_probs)
        cat("AUC (GAM fitted) =", ifelse(is.na(auc_val), "NA", round(auc_val, 4)), "\n")
        
        youden_cut <- find_youden_cut(model_data$outcome, model_data$Alb_by_NLR)
        cat("Breakpoint (Youden) =", ifelse(is.na(youden_cut), "NA", round(youden_cut, 3)), "\n")
        
        grid <- data.frame(
            Alb_by_NLR = seq(min(model_data$Alb_by_NLR, na.rm = TRUE),
                             max(model_data$Alb_by_NLR, na.rm = TRUE),
                             length.out = 200)
        )
        pr <- predict(gam_model, newdata = grid, type = "link", se.fit = TRUE)
        grid <- grid %>%
            mutate(
                prob_fit   = plogis(pr$fit),
                prob_lower = plogis(pr$fit - 1.96 * pr$se.fit),
                prob_upper = plogis(pr$fit + 1.96 * pr$se.fit)
            )
        
        bp <- NULL
        if (!is.na(youden_cut)) {
            bp <- grid %>%
                mutate(d = abs(Alb_by_NLR - youden_cut)) %>%
                filter(d == min(d, na.rm = TRUE)) %>%
                slice(1)
        }
        
        p <- ggplot(grid, aes(x = Alb_by_NLR)) +
            geom_line(aes(y = prob_lower), color = "black", linetype = "dotted") +
            geom_line(aes(y = prob_upper), color = "black", linetype = "dotted") +
            geom_line(aes(y = prob_fit), color = "red", linewidth = 1.1) +
            geom_rug(data = filter(model_data, outcome == 0),
                     inherit.aes = FALSE, aes(x = Alb_by_NLR),
                     sides = "b", alpha = 0.4, color = "steelblue") +
            geom_rug(data = filter(model_data, outcome == 1),
                     inherit.aes = FALSE, aes(x = Alb_by_NLR),
                     sides = "t", alpha = 0.4, color = "firebrick") +
            labs(title = panel_title, x = "Albumin/NLR", y = y_label) +
            theme_classic(base_size = 14) +
            theme(
                plot.title = element_text(face = "bold", size = 12),
                axis.title = element_text(face = "bold"),
                axis.text  = element_text(face = "bold", color = "black"),
                axis.line  = element_line(color = "black", linewidth = 0.5)
            ) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
        
        if (!is.null(bp) && nrow(bp) == 1) {
            p <- p +
                geom_vline(xintercept = bp$Alb_by_NLR, linetype = "dashed", color = "black") +
                geom_point(data = bp, aes(y = prob_fit), color = "black", size = 3.5, shape = 18) +
                geom_text(
                    data = bp,
                    aes(y = prob_fit, label = paste0("Youden = ", round(bp$Alb_by_NLR, 3))),
                    vjust = -1.2, fontface = "bold", color = "black"
                )
        }
        
        p
    }
    
    # ==============================================================================
    # MAIN (NO train/validation)
    # ==============================================================================
    stopifnot(exists("data"))
    
    data_cleaned  <- apply_exclusion_criteria(data)
    data_featured <- create_features(data_cleaned)
    
    # A) MASLD
    df_masld <- data_featured %>%
        filter(!is.na(masld_binary), masld_binary %in% c(0, 1),
               !is.na(Alb_by_NLR), is.finite(Alb_by_NLR))
    
    # B) SLF within MASLD
    df_slf <- data_featured %>%
        filter(masld_group == "MASLD",
               !is.na(slf_binary), slf_binary %in% c(0, 1),
               !is.na(Alb_by_NLR), is.finite(Alb_by_NLR))
    
    cat("--- PHÂN TÍCH GAM TRÊN TOÀN BỘ DỮ LIỆU (KHÔNG CHIA TRAIN/VALIDATION) ---\n")
    cat("MASLD N=", nrow(df_masld), "\n", sep = "")
    cat("SLF (within MASLD) N=", nrow(df_slf), "\n\n", sep = "")
    
    pA <- run_gam_panel(df_masld, "masld_binary", "A: MASLD (Full eligible cohort)", "Probability of MASLD")
    pB <- run_gam_panel(df_slf,   "slf_binary",   "B: SLF (MASLD cohort only)",      "Probability of SLF")
    
    combined <- pA | pB
    print(combined)
    
    cat("\n--- HOÀN TẤT ---\n")
    
    ```
    
- **Two picewise regression model (tự chạy điểm K luôn)**
    
    ```r
    # ==============================================================================
    # END-TO-END: Threshold (piecewise) logistic regression + LRT + deltaMethod + GT table
    # - NO train/validation split
    # - NO CURI, NO hs-CRP, NO uric acid, NO WBC
    # - Predictor: Alb_by_NLR = Albumin(g/dL) / NLR
    # - NLR priority: counts (LBDNENO/LBDLYMNO). Fallback: % (LBXNEPCT/LBXLYPCT) if counts missing
    # - Albumin priority: LBXSAL (g/dL). Fallback: LBDSALSI/10 (g/L -> g/dL)
    # - Exclusion key measurements: LUXCAPM, LUXSMED, Albumin + (Neut/Lymph counts OR %)
    # - Uses your required columns:
    #     LBXLYPCT  % lymphocyte
    #     LBXNEPCT  % neutrophil
    #     LBDLYMNO  lymphocyte (absolute)
    #     LBDNENO   neutrophil (absolute)
    #     LBXSAL    Albumin (g/dL)
    #     LBDSALSI  Albumin (g/L)
    # - Outcomes:
    #     MASLD: CAP>=263 + >=1 cardiometabolic risk (+ light drinker logic preserved)
    #     SLF: within MASLD only, LSM>=8.0
    # - K breakpoints are estimated from ROC-Youden on Alb_by_NLR (no manual placeholders)
    # ==============================================================================
    
    suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(tidyverse)
      library(gt)
      library(broom)
      library(lmtest)
      library(car)
    })
    
    # ----------------------------
    # Helpers
    # ----------------------------
    ensure_vars <- function(df, vars) {
      miss <- setdiff(vars, names(df))
      if (length(miss) > 0) for (v in miss) df[[v]] <- NA
      df
    }
    
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    format_p_value <- function(p) {
      if (is.na(p)) return(NA_character_)
      ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    }
    
    # =========================
    # 1) Alcohol consumption
    # =========================
    calculate_alcohol_consumption <- function(data) {
      data <- ensure_vars(data, c("ALQ121", "ALQ130"))
    
      data %>%
        mutate(
          ALQ121 = safe_num(ALQ121),
          ALQ130 = safe_num(ALQ130),
          days_per_week_alcohol = case_when(
            is.na(ALQ121) ~ NA_real_,
            ALQ121 == 0 ~ 0,
            ALQ121 == 1 ~ 7,
            ALQ121 == 2 ~ 6,
            ALQ121 == 3 ~ 3.5,
            ALQ121 == 4 ~ 2,
            ALQ121 == 5 ~ 1,
            ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
            ALQ121 == 7 ~ 1 / (30.4375 / 7),
            ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
            ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
            ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
            TRUE ~ NA_real_
          ),
          ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
          weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
        )
    }
    
    # =========================
    # 2) Exclusion criteria (key measurements = Alb/NLR inputs + FibroScan)
    # =========================
    apply_exclusion_criteria <- function(data) {
      cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ ---\n")
      initial_rows_total <- nrow(data)
    
      # drop fully-NA cols
      data <- data %>% dplyr::select(where(~ any(!is.na(.x))))
    
      # ensure columns (NO hs-CRP, NO uric acid)
      req <- c(
        "RIAGENDR","RIDAGEYR","LUXSMED","LUXCAPM","ALQ121","ALQ130",
        "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"
      )
      data <- ensure_vars(data, req)
    
      # 1) Viral hepatitis
      cat("1. Loại bỏ các trường hợp viêm gan virus...\n")
      hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR", "LBDHCV")
      hep_vars_exist <- intersect(hep_vars, names(data))
    
      rows_before <- nrow(data)
      if (length(hep_vars_exist) > 0) {
        data <- data %>%
          mutate(across(all_of(hep_vars_exist), safe_num)) %>%
          filter(rowSums(across(all_of(hep_vars_exist), ~ (.x == 1)), na.rm = TRUE) == 0)
      }
      cat("-> Đã loại bỏ", rows_before - nrow(data), "người do có bằng chứng viêm gan virus.\n\n")
    
      # 2) Drop NA on key measurements
      cat("2. Loại bỏ thiếu dữ liệu ở các biến đo lường chính (FibroScan + Alb/NLR inputs)...\n")
      rows_before <- nrow(data)
    
      data <- data %>%
        mutate(
          RIAGENDR = safe_num(RIAGENDR),
          RIDAGEYR = safe_num(RIDAGEYR),
          LUXSMED  = safe_num(LUXSMED),
          LUXCAPM  = safe_num(LUXCAPM),
    
          LBXSAL   = safe_num(LBXSAL),
          LBDSALSI = safe_num(LBDSALSI),
    
          LBDNENO  = safe_num(LBDNENO),
          LBDLYMNO = safe_num(LBDLYMNO),
          LBXNEPCT = safe_num(LBXNEPCT),
          LBXLYPCT = safe_num(LBXLYPCT),
    
          Albumin_g_dL_tmp = case_when(
            !is.na(LBXSAL) ~ LBXSAL,
            is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
            TRUE ~ NA_real_
          ),
          has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO) & LBDLYMNO > 0,
          has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT) & LBXLYPCT > 0
        ) %>%
        drop_na(LUXSMED, LUXCAPM) %>%
        filter(!is.na(Albumin_g_dL_tmp)) %>%
        filter(has_counts | has_pct)
    
      cat("-> Đã loại bỏ", rows_before - nrow(data), "hàng do thiếu FibroScan/Albumin/Neut-Lymph.\n\n")
    
      # 3) Alcohol
      cat("3. Loại bỏ uống rượu vượt ngưỡng...\n")
      data <- calculate_alcohol_consumption(data)
      rows_before <- nrow(data)
      data <- data %>%
        filter(
          is.na(weekly_alcohol_grams) |
            (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
            (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
        )
      cat("-> Đã loại bỏ", rows_before - nrow(data), "hàng do rượu vượt ngưỡng.\n\n")
    
      # 4) Age >= 18
      cat("4. Loại bỏ <18 tuổi...\n")
      rows_before <- nrow(data)
      data <- data %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
      cat("-> Đã loại bỏ", rows_before - nrow(data), "hàng do tuổi < 18.\n\n")
    
      cat("--- HOÀN TẤT LOẠI TRỪ. Tổng loại:", initial_rows_total - nrow(data),
          "Còn lại:", nrow(data), "---\n\n")
      data
    }
    
    # =========================
    # 3) Feature engineering (MASLD/SLF + Alb/NLR + covariates)
    # =========================
    create_features <- function(data) {
      cat("--- BẮT ĐẦU: TẠO CÁC BIẾN PHÂN TÍCH ---\n")
    
      cols_to_check <- c(
        "LBXSGL","SMQ020","SMQ040","DIQ010","DIQ050","DIQ070","BPQ040A",
        "MCQ160C","MCQ160E","MCQ160D","KIQ022","KIQ025","MCQ220","BPQ020",
        "LBDSCRSI","LBXSCR","DMDMARTZ","DMDEDUC2","RIDRETH1",
        "BMXBMI","BMXWAIST","LBXGLU","LBXGH",
        "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
        "LBDTRSI","BPQ090D","LBDHDDSI",
        "LUXCAPM","LUXSMED",
        # required Alb/NLR inputs
        "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI",
        "weekly_alcohol_grams",
        "RIAGENDR","RIDAGEYR"
      )
      data <- ensure_vars(data, cols_to_check)
    
      d <- data %>%
        mutate(across(all_of(cols_to_check), safe_num)) %>%
        mutate(
          # ---------------- MASLD risks ----------------
          has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
    
          risk_bmi_waist = if_else(
            is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
            NA,
            (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)
          ),
    
          risk_glucose_diabetes = if_else(
            is.na(LBXGLU) & is.na(LBXGH) &
              (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
              (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
              (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
            NA,
            (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
              (if_else(is.na(LBXGH),  FALSE, LBXGH  >= 5.7)) |
              (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
              (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
              (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
          ),
    
          any_bp_measurement_high =
            (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
            (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
            (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
    
          risk_blood_pressure = if_else(
            is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
            NA,
            any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
          ),
    
          risk_triglycerides = if_else(
            is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
            NA,
            (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          ),
    
          risk_low_hdl = if_else(
            is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))),
            NA,
            ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
               (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          ),
    
          num_cardiometabolic_risks = rowSums(
            cbind(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
            na.rm = TRUE
          ),
    
          has_one_plus_cardiometabolic_risk = if_else(
            is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
              is.na(risk_triglycerides) & is.na(risk_low_hdl),
            NA, num_cardiometabolic_risks > 0
          ),
    
          is_light_drinker_final = case_when(
            is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
            (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
            !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
            TRUE ~ NA
          ),
    
          masld_group = case_when(
            is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
            has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
            has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
            has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
            TRUE ~ NA_character_
          ),
          masld_binary = if_else(masld_group == "MASLD", 1L, 0L),
    
          slf_group = case_when(
            masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ "SLF",
            masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED <  8.0 ~ "non-SLF",
            TRUE ~ NA_character_
          ),
          slf_binary = case_when(
            masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ 1L,
            masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED <  8.0 ~ 0L,
            TRUE ~ NA_integer_
          ),
    
          # ---------------- Predictor Alb/NLR ----------------
          Albumin_g_dL = case_when(
            !is.na(LBXSAL) ~ LBXSAL,
            is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
            TRUE ~ NA_real_
          ),
          Albumin_g_L = case_when(
            !is.na(LBDSALSI) ~ LBDSALSI,
            is.na(LBDSALSI) & !is.na(LBXSAL) ~ LBXSAL * 10,
            TRUE ~ NA_real_
          ),
    
          Neutrophil_abs = LBDNENO,
          Lymphocyte_abs = LBDLYMNO,
          Neutrophil_pct = LBXNEPCT,
          Lymphocyte_pct = LBXLYPCT,
    
          NLR = case_when(
            !is.na(Neutrophil_abs) & !is.na(Lymphocyte_abs) & Lymphocyte_abs > 0 ~ Neutrophil_abs / Lymphocyte_abs,
            (is.na(Neutrophil_abs) | is.na(Lymphocyte_abs)) &
              !is.na(Neutrophil_pct) & !is.na(Lymphocyte_pct) & Lymphocyte_pct > 0 ~ Neutrophil_pct / Lymphocyte_pct,
            TRUE ~ NA_real_
          ),
    
          Alb_by_NLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                               Albumin_g_dL / NLR, NA_real_)
        ) %>%
        mutate(
          # eGFR for CKD covariate
          eGFR = {
            creatinine_mg_dl <- case_when(
              "LBXSCR" %in% names(.) ~ safe_num(LBXSCR),
              "LBDSCRSI" %in% names(.) ~ safe_num(LBDSCRSI) / 88.4,
              TRUE ~ NA_real_
            )
            if_else(
              !is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
              175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                if_else(RIAGENDR == 2, 0.742, 1) *
                if_else(RIDRETH1 == 4, 1.212, 1),
              NA_real_
            )
          },
    
          age = RIDAGEYR,
          gender = factor(RIAGENDR, labels = c("Male", "Female")),
          race = factor(case_when(
            RIDRETH1 == 1 ~ "Mexican American",
            RIDRETH1 == 2 ~ "Other Hispanic",
            RIDRETH1 == 3 ~ "Non-Hispanic White",
            RIDRETH1 == 4 ~ "Non-Hispanic Black",
            RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
            TRUE ~ "Other Race/Multi-Racial"
          )),
          marital_status = factor(case_when(
            DMDMARTZ %in% c(1, 6) ~ "Married/Living with partner",
            DMDMARTZ %in% c(2, 3, 4) ~ "Widowed/Divorced/Separated",
            DMDMARTZ == 5 ~ "Never married",
            TRUE ~ NA_character_
          )),
          education = factor(case_when(
            DMDEDUC2 == 1 ~ "Less than 9th grade",
            DMDEDUC2 == 2 ~ "9-11th grade",
            DMDEDUC2 == 3 ~ "High school graduate",
            DMDEDUC2 == 4 ~ "Some college or AA degree",
            DMDEDUC2 == 5 ~ "College graduate or above",
            TRUE ~ NA_character_
          )),
          smoking_status = factor(case_when(
            SMQ020 == 2 ~ "Non-smoker",
            SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
            SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
            TRUE ~ NA_character_
          ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
          obesity = factor(if_else(BMXBMI >= 30, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
          diabetes_mellitus = factor(if_else(
            DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
              (!is.na(LBXGH) & LBXGH >= 6.5) |
              (!is.na(LBXSGL) & LBXSGL >= 126),
            "Yes", "No", missing = "No"
          ), levels = c("No", "Yes")),
          hypertension = factor(if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
          chronic_kidney_disease = factor(if_else((KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60)),
                                                  "Yes", "No", missing = "No"), levels = c("No", "Yes")),
          cancer = factor(if_else(MCQ220 == 1, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
          cardiovascular_disease = factor(if_else(MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1,
                                                 "Yes", "No", missing = "No"), levels = c("No", "Yes"))
        )
    
      cat("--- HOÀN TẤT TẠO BIẾN. ---\n\n")
      d
    }
    
    # =========================
    # 4) Youden breakpoint on Alb/NLR
    # =========================
    youden_cutoff <- function(df, outcome, predictor) {
      ok <- !is.na(df[[outcome]]) & df[[outcome]] %in% c(0, 1) & !is.na(df[[predictor]]) & is.finite(df[[predictor]])
      d <- df[ok, , drop = FALSE]
      if (nrow(d) < 50) return(NA_real_)
      if (n_distinct(d[[outcome]]) < 2) return(NA_real_)
    
      roc_obj <- tryCatch(
        pROC::roc(response = d[[outcome]], predictor = d[[predictor]], levels = c(0, 1), direction = "<", quiet = TRUE),
        error = function(e) NULL
      )
      if (is.null(roc_obj)) return(NA_real_)
    
      thr <- tryCatch(
        as.numeric(pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold")[1]),
        error = function(e) NA_real_
      )
      thr
    }
    
    # =========================
    # 5) Piecewise logistic regression + LRT + deltaMethod
    # =========================
    perform_piecewise_analysis <- function(dataset, k_masld, k_slf, scaling_factor = 1) {
    
      results <- list()
      results$scaling_factor <- scaling_factor
    
      covs_all <- c(
        "age", "gender", "race", "marital_status", "education",
        "smoking_status", "obesity", "diabetes_mellitus", "hypertension",
        "chronic_kidney_disease", "cancer", "cardiovascular_disease"
      )
      covs <- covs_all %>%
        keep(~ .x %in% names(dataset) && n_distinct(dataset[[.x]], na.rm = TRUE) >= 2)
    
      cat("--- Covariates used:", paste(covs, collapse = ", "), "\n\n")
    
      k_masld_scaled <- k_masld / scaling_factor
      k_slf_scaled   <- k_slf   / scaling_factor
    
      # ---- MASLD ----
      dm <- dataset %>%
        select(masld_binary, Alb_by_NLR, all_of(covs)) %>%
        filter(!is.na(masld_binary), masld_binary %in% c(0, 1),
               !is.na(Alb_by_NLR), is.finite(Alb_by_NLR)) %>%
        na.omit() %>%
        mutate(Alb_by_NLR_scaled = Alb_by_NLR / scaling_factor)
    
      results$n_masld <- nrow(dm)
    
      if (nrow(dm) > length(covs) + 10 && is.finite(k_masld_scaled)) {
        dm <- dm %>%
          mutate(
            segment1 = Alb_by_NLR_scaled,
            segment2 = pmax(0, Alb_by_NLR_scaled - k_masld_scaled)
          )
    
        f_pw <- as.formula(paste("masld_binary ~ segment1 + segment2 +", paste(covs, collapse = " + ")))
        f_ln <- as.formula(paste("masld_binary ~ Alb_by_NLR_scaled +", paste(covs, collapse = " + ")))
    
        pw <- glm(f_pw, data = dm, family = binomial())
        ln <- glm(f_ln, data = dm, family = binomial())
    
        s <- summary(pw)
    
        results$or1_masld <- exp(coef(pw)["segment1"])
        results$ci1_masld <- exp(confint.default(pw)["segment1", ])
        results$p1_masld  <- s$coefficients["segment1", "Pr(>|z|)"]
    
        # slope after K: exp(segment1 + segment2)
        dm2 <- car::deltaMethod(pw, "exp(segment1 + segment2)")
        results$or2_masld <- unname(dm2$Estimate)
        results$ci2_masld <- c(unname(dm2$`2.5 %`), unname(dm2$`97.5 %`))
        # p for change in slope (segment2)
        results$p2_masld  <- s$coefficients["segment2", "Pr(>|z|)"]
    
        lrt <- lmtest::lrtest(ln, pw)
        results$lrt_masld <- lrt$`Pr(>Chisq)`[2]
      } else {
        results$or1_masld <- results$ci1_masld <- results$p1_masld <- NA
        results$or2_masld <- results$ci2_masld <- results$p2_masld <- NA
        results$lrt_masld <- NA
      }
    
      # ---- SLF (within MASLD) ----
      ds <- dataset %>%
        filter(masld_group == "MASLD") %>%
        select(slf_binary, Alb_by_NLR, all_of(covs)) %>%
        filter(!is.na(slf_binary), slf_binary %in% c(0, 1),
               !is.na(Alb_by_NLR), is.finite(Alb_by_NLR)) %>%
        na.omit() %>%
        mutate(Alb_by_NLR_scaled = Alb_by_NLR / scaling_factor)
    
      results$n_slf <- nrow(ds)
    
      if (nrow(ds) > length(covs) + 10 && is.finite(k_slf_scaled)) {
        ds <- ds %>%
          mutate(
            segment1 = Alb_by_NLR_scaled,
            segment2 = pmax(0, Alb_by_NLR_scaled - k_slf_scaled)
          )
    
        f_pw <- as.formula(paste("slf_binary ~ segment1 + segment2 +", paste(covs, collapse = " + ")))
        f_ln <- as.formula(paste("slf_binary ~ Alb_by_NLR_scaled +", paste(covs, collapse = " + ")))
    
        pw <- glm(f_pw, data = ds, family = binomial())
        ln <- glm(f_ln, data = ds, family = binomial())
    
        s <- summary(pw)
    
        results$or1_slf <- exp(coef(pw)["segment1"])
        results$ci1_slf <- exp(confint.default(pw)["segment1", ])
        results$p1_slf  <- s$coefficients["segment1", "Pr(>|z|)"]
    
        dm2 <- car::deltaMethod(pw, "exp(segment1 + segment2)")
        results$or2_slf <- unname(dm2$Estimate)
        results$ci2_slf <- c(unname(dm2$`2.5 %`), unname(dm2$`97.5 %`))
        results$p2_slf  <- s$coefficients["segment2", "Pr(>|z|)"]
    
        lrt <- lmtest::lrtest(ln, pw)
        results$lrt_slf <- lrt$`Pr(>Chisq)`[2]
      } else {
        results$or1_slf <- results$ci1_slf <- results$p1_slf <- NA
        results$or2_slf <- results$ci2_slf <- results$p2_slf <- NA
        results$lrt_slf <- NA
      }
    
      results
    }
    
    # =========================
    # 6) GT results table
    # =========================
    create_results_table <- function(results, k_masld, k_slf, table_title) {
    
      or_ci <- function(or, ci) {
        if (!is.null(or) && !is.null(ci) && all(is.finite(c(or, ci)))) {
          sprintf("%.2f (%.2f, %.2f)", or, ci[1], ci[2])
        } else {
          "N/A"
        }
      }
    
      df <- tibble(
        Models = c(
          "N",
          "K (Youden)",
          paste0("Alb/NLR < K (per ", results$scaling_factor, " unit)"),
          paste0("Alb/NLR > K (per ", results$scaling_factor, " unit)"),
          "Likelihood ratio test (piecewise vs linear)"
        ),
        `OR (95% CI)_masld` = c(
          results$n_masld,
          ifelse(is.na(k_masld), "NA", sprintf("%.3f", k_masld)),
          or_ci(results$or1_masld, results$ci1_masld),
          or_ci(results$or2_masld, results$ci2_masld),
          ""
        ),
        `P-value_masld` = c(
          "", "", format_p_value(results$p1_masld), format_p_value(results$p2_masld), format_p_value(results$lrt_masld)
        ),
        `OR (95% CI)_slf` = c(
          results$n_slf,
          ifelse(is.na(k_slf), "NA", sprintf("%.3f", k_slf)),
          or_ci(results$or1_slf, results$ci1_slf),
          or_ci(results$or2_slf, results$ci2_slf),
          ""
        ),
        `P-value_slf` = c(
          "", "", format_p_value(results$p1_slf), format_p_value(results$p2_slf), format_p_value(results$lrt_slf)
        )
      )
    
      df %>%
        gt() %>%
        tab_header(title = table_title) %>%
        tab_spanner(label = "MASLD", columns = c(`OR (95% CI)_masld`, `P-value_masld`)) %>%
        tab_spanner(label = "Significant Clinical Fibrosis", columns = c(`OR (95% CI)_slf`, `P-value_slf`)) %>%
        cols_label(
          Models = "",
          `OR (95% CI)_masld` = "OR (95% CI)", `P-value_masld` = "P-value",
          `OR (95% CI)_slf`   = "OR (95% CI)", `P-value_slf`   = "P-value"
        ) %>%
        cols_align(align = "left", columns = Models) %>%
        cols_align(align = "center", columns = -Models) %>%
        tab_style(
          style = cell_text(weight = "bold"),
          locations = cells_body(
            columns = c(`P-value_masld`, `P-value_slf`),
            rows = grepl("^<0\\.001$", `P-value_masld`) |
              suppressWarnings(as.numeric(`P-value_masld`)) < 0.05 |
              grepl("^<0\\.001$", `P-value_slf`) |
              suppressWarnings(as.numeric(`P-value_slf`)) < 0.05
          )
        ) %>%
        tab_source_note(md(
          "Adjusted for age, gender, race, marital status, education, smoking status, obesity, diabetes mellitus, hypertension, chronic kidney disease, cancer, and cardiovascular disease.<br>
          Predictor Alb/NLR computed from: LBDNENO/LBDLYMNO (primary) or LBXNEPCT/LBXLYPCT (fallback); Albumin from LBXSAL (primary) or LBDSALSI/10 (fallback)."
        ))
    }
    
    # ==============================================================================
    # MAIN
    # ==============================================================================
    stopifnot(exists("MASLD") || exists("data"))
    if (!exists("data")) data <- MASLD
    
    data_cleaned  <- apply_exclusion_criteria(data)
    data_featured <- create_features(data_cleaned)
    
    # Estimate K from ROC-Youden on Alb_by_NLR (NOT CURI)
    K_MASLD <- youden_cutoff(data_featured, "masld_binary", "Alb_by_NLR")
    K_SLF   <- youden_cutoff(dplyr::filter(data_featured, masld_group == "MASLD"), "slf_binary", "Alb_by_NLR")
    
    cat("--- Youden breakpoints (Alb/NLR) ---\n")
    cat("K_MASLD =", ifelse(is.na(K_MASLD), "NA", round(K_MASLD, 3)), "\n")
    cat("K_SLF   =", ifelse(is.na(K_SLF),   "NA", round(K_SLF, 3)), "\n\n")
    
    cat("--- Threshold (piecewise) analysis trên TOÀN BỘ DỮ LIỆU ---\n")
    full_results <- perform_piecewise_analysis(
      dataset = data_featured,
      k_masld = K_MASLD,
      k_slf   = K_SLF,
      scaling_factor = 1
    )
    
    final_table <- create_results_table(
      results = full_results,
      k_masld = K_MASLD,
      k_slf   = K_SLF,
      table_title = "Table. Threshold Effect Analysis (Full cohort)"
    )
    
    print(final_table)
    
    cat("\n--- HOÀN TẤT ---\n")
    
    ```
    
- **Two-picewise regression model train**
    
    ```r
    # ==============================================================================
    # END-TO-END: Threshold (piecewise) logistic regression + LRT + deltaMethod + GT table
    # - NO train/validation
    # - Predictor: Alb_by_NLR = Albumin(g/dL) / NLR ; NLR = LBDNENO / LBDLYMNO
    # - REQUIRED cols (no WBC):
    #     LBXLYPCT, LBXNEPCT, LBDLYMNO, LBDNENO, LBXSAL or LBDSALSI
    # - K points:
    #     K_MASLD = 0.606
    #     K_SLF   = 5.642
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(tidyverse)
        library(gt)
        library(broom)
        library(lmtest)
        library(car)   # deltaMethod
    })
    
    # ----------------------------
    # Helper: ensure vars exist (create NA cols if missing)
    # ----------------------------
    ensure_vars <- function(df, vars) {
        miss <- setdiff(vars, names(df))
        if (length(miss) > 0) for (v in miss) df[[v]] <- NA
        df
    }
    
    # ----------------------------
    # Helper: pick albumin (g/dL) from LBXSAL or LBDSALSI
    # ----------------------------
    get_albumin_gdl <- function(df) {
        alb_gdl <- rep(NA_real_, nrow(df))
        if ("LBXSAL" %in% names(df)) {
            alb_gdl <- suppressWarnings(as.numeric(df$LBXSAL))
        }
        if ("LBDSALSI" %in% names(df)) {
            alb_si <- suppressWarnings(as.numeric(df$LBDSALSI))      # g/L
            alb_si_gdl <- alb_si / 10                                # -> g/dL
            alb_gdl <- ifelse(is.na(alb_gdl) & !is.na(alb_si_gdl), alb_si_gdl, alb_gdl)
        }
        alb_gdl
    }
    
    # =========================
    # 1) Alcohol consumption
    # =========================
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_vars(data, c("ALQ121", "ALQ130"))
        
        data %>%
            mutate(
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, as.numeric(ALQ130)),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # =========================
    # 2) Exclusion criteria (NO hsCRP, NO uric acid)
    # =========================
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ ---\n")
        initial_rows_total <- nrow(data)
        
        # remove all-NA columns
        data <- data[, colSums(!is.na(data)) > 0]
        
        # ensure vars used in exclusion exist
        data <- ensure_vars(
            data,
            c("RIAGENDR","RIDAGEYR","LUXSMED","LUXCAPM",
              "ALQ121","ALQ130",
              "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO",
              "LBXSAL","LBDSALSI")
        )
        
        # 1) Viral hepatitis
        cat("1. Loại bỏ viêm gan virus...\n")
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR")
        hep_vars_exist <- intersect(hep_vars, names(data))
        rows_before <- nrow(data)
        if (length(hep_vars_exist) > 0) {
            data <- data %>% filter(rowSums(across(all_of(hep_vars_exist), ~ (.x == 1)), na.rm = TRUE) == 0)
        }
        cat("-> Loại:", rows_before - nrow(data), "\n\n")
        
        # 2) Drop NA on key measurements (FibroScan + Alb/NLR inputs)
        cat("2. Loại thiếu dữ liệu biến chính (FibroScan + Alb/NLR)...\n")
        rows_before <- nrow(data)
        
        data$Albumin_g_dl_tmp <- get_albumin_gdl(data)
        
        vars_must <- c("LUXSMED","LUXCAPM","LBDNENO","LBDLYMNO","Albumin_g_dl_tmp")
        data <- data %>% drop_na(all_of(vars_must))
        
        cat("-> Loại:", rows_before - nrow(data), "\n\n")
        data$Albumin_g_dl_tmp <- NULL
        
        # 3) Alcohol threshold
        cat("3. Loại uống rượu vượt ngưỡng...\n")
        data <- calculate_alcohol_consumption(data)
        rows_before <- nrow(data)
        data <- data %>%
            filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        cat("-> Loại:", rows_before - nrow(data), "\n\n")
        
        # 4) Age >= 18
        cat("4. Loại tuổi < 18...\n")
        rows_before <- nrow(data)
        data <- data %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        cat("-> Loại:", rows_before - nrow(data), "\n\n")
        
        cat("--- HOÀN TẤT LOẠI TRỪ. Tổng loại:", initial_rows_total - nrow(data),
            " | Còn lại:", nrow(data), "---\n\n")
        data
    }
    
    # =========================
    # 3) Feature engineering (MASLD/SLF + Alb/NLR + covariates)
    # =========================
    create_features <- function(data) {
        cat("--- BẮT ĐẦU: TẠO CÁC BIẾN PHÂN TÍCH ---\n")
        
        cols_to_check <- c(
            "RIAGENDR","RIDAGEYR","RIDRETH1","DMDMARTZ","DMDEDUC2",
            "SMQ020","SMQ040",
            "BMXBMI","BMXWAIST",
            "LBXGLU","LBXGH","LBXSGL",
            "DIQ010","DIQ050","DIQ070",
            "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
            "BPQ040A","BPQ090D","BPQ020",
            "LBDTRSI","LBDHDDSI",
            "MCQ160C","MCQ160E","MCQ160D","MCQ220",
            "KIQ022","KIQ025",
            "LBDSCRSI","LBXSCR",
            "LUXCAPM","LUXSMED",
            "weekly_alcohol_grams",
            "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO",
            "LBXSAL","LBDSALSI"
        )
        data <- ensure_vars(data, cols_to_check)
        
        data <- data %>%
            mutate(
                Albumin_g_dl   = get_albumin_gdl(.),
                Neutrophil_abs = suppressWarnings(as.numeric(LBDNENO)),
                Lymphocyte_abs = suppressWarnings(as.numeric(LBDLYMNO)),
                NLR = if_else(!is.na(Neutrophil_abs) & !is.na(Lymphocyte_abs) & Lymphocyte_abs > 0,
                              Neutrophil_abs / Lymphocyte_abs, NA_real_),
                Alb_by_NLR = if_else(!is.na(Albumin_g_dl) & !is.na(NLR) & NLR > 0,
                                     Albumin_g_dl / NLR, NA_real_)
            )
        
        data_processed <- data %>%
            mutate(
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                
                risk_bmi_waist = if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                    NA,
                    (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)
                ),
                
                risk_glucose_diabetes = if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (if_else(is.na(LBXGH),  FALSE, LBXGH  >= 5.7)) |
                        (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                
                any_bp_measurement_high =
                    (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                
                risk_blood_pressure = if_else(
                    is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
                    NA,
                    any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                
                risk_triglycerides = if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
                    NA,
                    (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                risk_low_hdl = if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))),
                    NA,
                    ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                )
            ) %>%
            mutate(
                num_cardiometabolic_risks = rowSums(
                    cbind(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
                    na.rm = TRUE
                ),
                has_one_plus_cardiometabolic_risk = if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA, num_cardiometabolic_risks > 0
                ),
                is_light_drinker_final = case_when(
                    is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
                    (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                masld_group = case_when(
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
                    TRUE ~ NA_character_
                ),
                masld_binary = if_else(masld_group == "MASLD", 1L, 0L),
                slf_group = case_when(
                    masld_group == "MASLD" & LUXSMED >= 8.0 ~ "SLF",
                    masld_group == "MASLD" & LUXSMED < 8.0  ~ "non-SLF",
                    TRUE ~ NA_character_
                ),
                slf_binary = case_when(
                    masld_group == "MASLD" & LUXSMED >= 8.0 ~ 1L,
                    masld_group == "MASLD" & LUXSMED < 8.0  ~ 0L,
                    TRUE ~ NA_integer_
                )
            ) %>%
            mutate(
                eGFR = {
                    creatinine_mg_dl <- case_when(
                        "LBXSCR" %in% names(.) ~ suppressWarnings(as.numeric(LBXSCR)),
                        "LBDSCRSI" %in% names(.) ~ suppressWarnings(as.numeric(LBDSCRSI)) / 88.4,
                        TRUE ~ NA_real_
                    )
                    if_else(
                        !is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                        175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                            if_else(RIAGENDR == 2, 0.742, 1) *
                            if_else(RIDRETH1 == 4, 1.212, 1),
                        NA_real_
                    )
                },
                age = RIDAGEYR,
                gender = factor(RIAGENDR, labels = c("Male", "Female")),
                race = factor(case_when(
                    RIDRETH1 == 1 ~ "Mexican American",
                    RIDRETH1 == 2 ~ "Other Hispanic",
                    RIDRETH1 == 3 ~ "Non-Hispanic White",
                    RIDRETH1 == 4 ~ "Non-Hispanic Black",
                    RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
                    TRUE ~ "Other Race/Multi-Racial"
                )),
                marital_status = factor(case_when(
                    DMDMARTZ %in% c(1, 6) ~ "Married/Living with partner",
                    DMDMARTZ %in% c(2, 3, 4) ~ "Widowed/Divorced/Separated",
                    DMDMARTZ == 5 ~ "Never married",
                    TRUE ~ NA_character_
                )),
                education = factor(case_when(
                    DMDEDUC2 == 1 ~ "Less than 9th grade",
                    DMDEDUC2 == 2 ~ "9-11th grade",
                    DMDEDUC2 == 3 ~ "High school graduate",
                    DMDEDUC2 == 4 ~ "Some college or AA degree",
                    DMDEDUC2 == 5 ~ "College graduate or above",
                    TRUE ~ NA_character_
                )),
                smoking_status = factor(case_when(
                    SMQ020 == 2 ~ "Non-smoker",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
                    TRUE ~ NA_character_
                ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
                obesity = factor(if_else(BMXBMI >= 30, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
                diabetes_mellitus = factor(if_else(
                    DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                        (!is.na(LBXGH) & LBXGH >= 6.5) |
                        (!is.na(LBXSGL) & LBXSGL >= 126),
                    "Yes", "No", missing = "No"
                ), levels = c("No", "Yes")),
                hypertension = factor(if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
                chronic_kidney_disease = factor(if_else((KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60)),
                                                        "Yes", "No", missing = "No"), levels = c("No", "Yes")),
                cancer = factor(if_else(MCQ220 == 1, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
                cardiovascular_disease = factor(if_else(MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1,
                                                        "Yes", "No", missing = "No"), levels = c("No", "Yes"))
            )
        
        cat("--- HOÀN TẤT TẠO BIẾN. ---\n\n")
        data_processed
    }
    
    format_p_value <- function(p) {
        if (is.na(p)) return(NA_character_)
        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    }
    
    # =========================
    # Piecewise logistic regression
    # =========================
    perform_piecewise_analysis <- function(dataset, k_masld, k_slf, scaling_factor = 1) {
        
        results <- list()
        final_covariates <- c(
            "age", "gender", "race", "marital_status", "education",
            "smoking_status", "obesity", "diabetes_mellitus", "hypertension",
            "chronic_kidney_disease", "cancer", "cardiovascular_disease"
        )
        final_covariates <- final_covariates %>%
            keep(~ .x %in% names(dataset) && n_distinct(dataset[[.x]], na.rm = TRUE) >= 2)
        
        cat("--- Covariates used:", paste(final_covariates, collapse = ", "), "\n\n")
        
        results$scaling_factor <- scaling_factor
        k_masld_scaled <- k_masld / scaling_factor
        k_slf_scaled   <- k_slf   / scaling_factor
        
        # ---------- MASLD ----------
        model_data_masld <- dataset %>%
            dplyr::select(masld_binary, Alb_by_NLR, all_of(final_covariates)) %>%
            filter(!is.na(masld_binary), masld_binary %in% c(0, 1),
                   !is.na(Alb_by_NLR), is.finite(Alb_by_NLR)) %>%
            na.omit() %>%
            mutate(Alb_by_NLR_scaled = Alb_by_NLR / scaling_factor)
        
        results$n_masld <- nrow(model_data_masld)
        
        if (nrow(model_data_masld) > length(final_covariates) + 5) {
            
            model_data_masld <- model_data_masld %>%
                mutate(
                    segment1 = Alb_by_NLR_scaled,
                    segment2 = pmax(0, Alb_by_NLR_scaled - k_masld_scaled)
                )
            
            piecewise_formula_masld <- as.formula(paste(
                "masld_binary ~ segment1 + segment2 +",
                paste(final_covariates, collapse = " + ")
            ))
            
            linear_formula_masld <- as.formula(paste(
                "masld_binary ~ Alb_by_NLR_scaled +",
                paste(final_covariates, collapse = " + ")
            ))
            
            tryCatch({
                pw <- glm(piecewise_formula_masld, data = model_data_masld, family = binomial())
                lin <- glm(linear_formula_masld,   data = model_data_masld, family = binomial())
                
                summ <- summary(pw)
                
                results$or1_masld <- exp(coef(pw)["segment1"])
                results$ci1_masld <- exp(confint.default(pw)["segment1", ])
                results$p1_masld  <- summ$coefficients["segment1", "Pr(>|z|)"]
                
                or2_calc <- car::deltaMethod(pw, "exp(segment1 + segment2)")
                results$or2_masld <- unname(or2_calc$Estimate)
                results$ci2_masld <- c(unname(or2_calc$`2.5 %`), unname(or2_calc$`97.5 %`))
                results$p2_masld  <- summ$coefficients["segment2", "Pr(>|z|)"]
                
                lrt <- lmtest::lrtest(lin, pw)
                results$lrt_masld <- lrt$`Pr(>Chisq)`[2]
            }, error = function(e) cat("LỖI (MASLD):", conditionMessage(e), "\n"))
        }
        
        # ---------- SLF ----------
        model_data_slf <- dataset %>%
            filter(masld_group == "MASLD") %>%
            dplyr::select(slf_binary, Alb_by_NLR, all_of(final_covariates)) %>%
            filter(!is.na(slf_binary), slf_binary %in% c(0, 1),
                   !is.na(Alb_by_NLR), is.finite(Alb_by_NLR)) %>%
            na.omit() %>%
            mutate(Alb_by_NLR_scaled = Alb_by_NLR / scaling_factor)
        
        results$n_slf <- nrow(model_data_slf)
        
        if (nrow(model_data_slf) > length(final_covariates) + 5) {
            
            model_data_slf <- model_data_slf %>%
                mutate(
                    segment1 = Alb_by_NLR_scaled,
                    segment2 = pmax(0, Alb_by_NLR_scaled - k_slf_scaled)
                )
            
            piecewise_formula_slf <- as.formula(paste(
                "slf_binary ~ segment1 + segment2 +",
                paste(final_covariates, collapse = " + ")
            ))
            
            linear_formula_slf <- as.formula(paste(
                "slf_binary ~ Alb_by_NLR_scaled +",
                paste(final_covariates, collapse = " + ")
            ))
            
            tryCatch({
                pw <- glm(piecewise_formula_slf, data = model_data_slf, family = binomial())
                lin <- glm(linear_formula_slf,   data = model_data_slf, family = binomial())
                
                summ <- summary(pw)
                
                results$or1_slf <- exp(coef(pw)["segment1"])
                results$ci1_slf <- exp(confint.default(pw)["segment1", ])
                results$p1_slf  <- summ$coefficients["segment1", "Pr(>|z|)"]
                
                or2_calc <- car::deltaMethod(pw, "exp(segment1 + segment2)")
                results$or2_slf <- unname(or2_calc$Estimate)
                results$ci2_slf <- c(unname(or2_calc$`2.5 %`), unname(or2_calc$`97.5 %`))
                results$p2_slf  <- summ$coefficients["segment2", "Pr(>|z|)"]
                
                lrt <- lmtest::lrtest(lin, pw)
                results$lrt_slf <- lrt$`Pr(>Chisq)`[2]
            }, error = function(e) cat("LỖI (SLF):", conditionMessage(e), "\n"))
        }
        
        results
    }
    
    # =========================
    # GT table
    # =========================
    create_results_table <- function(results, k_masld, k_slf, table_title) {
        or_ci <- function(or, ci) {
            if (!is.null(or) && !is.null(ci) && all(is.finite(c(or, ci)))) {
                sprintf("%.2f (%.2f, %.2f)", or, ci[1], ci[2])
            } else {
                "N/A"
            }
        }
        
        results_df <- tibble(
            Models = c(
                "Số quan sát (N)",
                "Điểm gãy (K)",
                paste0("Albumin/NLR < K (per ", results$scaling_factor, " units)"),
                paste0("Albumin/NLR > K (per ", results$scaling_factor, " units)"),
                "Log-likelihood ratio test"
            ),
            `OR (95% CI)_masld` = c(
                results$n_masld,
                sprintf("%.3f", k_masld),
                or_ci(results$or1_masld, results$ci1_masld),
                or_ci(results$or2_masld, results$ci2_masld),
                ""
            ),
            `P-value_masld` = c(
                "", "", format_p_value(results$p1_masld), format_p_value(results$p2_masld), format_p_value(results$lrt_masld)
            ),
            `OR (95% CI)_slf` = c(
                results$n_slf,
                sprintf("%.3f", k_slf),
                or_ci(results$or1_slf, results$ci1_slf),
                or_ci(results$or2_slf, results$ci2_slf),
                ""
            ),
            `P-value_slf` = c(
                "", "", format_p_value(results$p1_slf), format_p_value(results$p2_slf), format_p_value(results$lrt_slf)
            )
        )
        
        results_df %>%
            gt() %>%
            tab_header(title = table_title) %>%
            tab_spanner(label = "MASLD", columns = c(`OR (95% CI)_masld`, `P-value_masld`)) %>%
            tab_spanner(label = "Significant Clinical Fibrosis", columns = c(`OR (95% CI)_slf`, `P-value_slf`)) %>%
            cols_label(
                Models = "Models",
                `OR (95% CI)_masld` = "OR (95% CI)", `P-value_masld` = "P-value",
                `OR (95% CI)_slf`   = "OR (95% CI)", `P-value_slf`   = "P-value"
            ) %>%
            cols_align(align = "center", columns = -Models) %>%
            cols_align(align = "left", columns = Models) %>%
            tab_source_note(md(
                "Adjusted for age, gender, race, marital status, education, smoking status, obesity, diabetes mellitus, hypertension, chronic kidney disease, cancer, and cardiovascular disease."
            ))
    }
    
    # ==============================================================================
    # MAIN
    # ==============================================================================
    stopifnot(exists("data") || exists("MASLD"))
    if (!exists("data")) data <- MASLD
    
    data_cleaned  <- apply_exclusion_criteria(data)
    data_featured <- create_features(data_cleaned)
    
    # K from your GAM+Youden output
    K_MASLD <- 0.606
    K_SLF   <- 5.642
    
    cat("--- THRESHOLD (PIECEWISE) ANALYSIS (FULL DATA) ---\n")
    cat("K_MASLD =", K_MASLD, " | K_SLF =", K_SLF, "\n\n")
    
    full_results <- perform_piecewise_analysis(
        dataset = data_featured,
        k_masld = K_MASLD,
        k_slf   = K_SLF,
        scaling_factor = 1
    )
    
    final_table <- create_results_table(
        results = full_results,
        k_masld = K_MASLD,
        k_slf   = K_SLF,
        table_title = "Table. Threshold Effect Analysis (Alb/NLR; full cohort)"
    )
    
    print(final_table)
    
    cat("\n--- HOÀN TẤT ---\n")
    
    ```
    
- **Subgroup**
    
    ```r
    # ====================== END-TO-END: 4-PANEL FOREST (NO N COLUMN, ALNR LABEL, CLEAR FORMAT) ======================
    # NO train/validation split
    # Exclusion key measurements: ONLY ALNR inputs (NO hs-CRP, NO uric acid, NO WBC)
    # Use exactly these columns for ALNR:
    #   LBXLYPCT  % lymphocyte
    #   LBXNEPCT  % neutrophil
    #   LBDLYMNO  lymphocyte (absolute)
    #   LBDNENO   neutrophil (absolute)
    #   LBXSAL    Albumin (g/dL)
    #   LBDSALSI  Albumin (g/L)
    # NLR priority: abs counts (LBDNENO/LBDLYMNO). Fallback: % (LBXNEPCT/LBXLYPCT) if counts missing
    # Albumin priority: LBXSAL, fallback: LBDSALSI/10
    # ==============================================================================================================
    
    suppressPackageStartupMessages({
        library(dplyr); library(tidyr); library(purrr); library(stringr)
        library(ggplot2); library(patchwork); library(scales)
        library(forcats); library(logistf)
    })
    
    # ----------------------------- 0) HELPERS -----------------------------
    ensure_cols <- function(df, cols){
        miss <- setdiff(cols, names(df))
        if (length(miss)) for (m in miss) df[[m]] <- NA
        df
    }
    
    num_safe <- function(x){
        suppressWarnings(as.numeric(x))
    }
    
    # Alb (g/dL): LBXSAL priority, else LBDSALSI (g/L) / 10
    get_albumin_gdl <- function(df){
        lbx <- num_safe(df$LBXSAL)
        lbs <- num_safe(df$LBDSALSI)
        ifelse(!is.na(lbx), lbx,
               ifelse(!is.na(lbs), lbs/10, NA_real_))
    }
    
    # NLR: priority abs counts; fallback % if abs missing
    get_nlr <- function(df){
        neu_abs <- num_safe(df$LBDNENO)
        lym_abs <- num_safe(df$LBDLYMNO)
        neu_pct <- num_safe(df$LBXNEPCT)
        lym_pct <- num_safe(df$LBXLYPCT)
        
        nlr_abs <- ifelse(!is.na(neu_abs) & !is.na(lym_abs) & lym_abs > 0, neu_abs/lym_abs, NA_real_)
        nlr_pct <- ifelse(!is.na(neu_pct) & !is.na(lym_pct) & lym_pct > 0, neu_pct/lym_pct, NA_real_)
        
        ifelse(!is.na(nlr_abs), nlr_abs, nlr_pct)
    }
    
    # ----------------------------- 1) EXCLUSION -----------------------------
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_cols(data, c("ALQ121","ALQ130"))
        data %>%
            mutate(
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, as.numeric(ALQ130)),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    apply_exclusion_criteria <- function(data) {
        
        # drop all-NA cols
        data <- data[, colSums(!is.na(data)) > 0, drop = FALSE]
        
        # ensure columns referenced anywhere below
        data <- ensure_cols(data, c(
            "RIAGENDR","RIDAGEYR",
            "LUXSMED","LUXCAPM",
            # ALNR inputs (EXACT set)
            "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI",
            "ALQ121","ALQ130"
        ))
        
        # viral hepatitis exclusion (if present)
        hep_vars <- c("RHD143","LBXHBC","LBDHBG","LBXHCR")
        hep_vars_exist <- intersect(hep_vars, names(data))
        if (length(hep_vars_exist) > 0) {
            data <- data %>% filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
        }
        
        # MUST-HAVE: FibroScan vars for defining MASLD/SLF + ALNR inputs available to compute
        data <- data %>%
            filter(!is.na(LUXSMED), !is.na(LUXCAPM))
        
        # Compute Alb + NLR using EXACT inputs, then require non-missing ALNR
        data <- data %>%
            mutate(
                Albumin_g_dl = get_albumin_gdl(cur_data()),
                NLR          = get_nlr(cur_data()),
                ALNR         = if_else(!is.na(Albumin_g_dl) & !is.na(NLR) & NLR > 0, Albumin_g_dl / NLR, NA_real_)
            ) %>%
            filter(!is.na(ALNR))
        
        # alcohol exclusion (skip if ALQ missing)
        data <- calculate_alcohol_consumption(data)
        if ("RIAGENDR" %in% names(data) && "weekly_alcohol_grams" %in% names(data)) {
            data <- data %>%
                filter(
                    is.na(weekly_alcohol_grams) |
                        (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                        (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
                )
        }
        
        # age >=18 (if present)
        if ("RIDAGEYR" %in% names(data)) data <- data %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        
        data
    }
    
    # ----------------------------- 2) FEATURES -----------------------------
    create_features <- function(data) {
        
        need <- c(
            "RIAGENDR","RIDAGEYR","RIDRETH3","RIDRETH1",
            "BMXBMI","BMXWAIST","LBXGLU","LBXGH","DIQ010","DIQ050","DIQ070",
            "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3","BPQ040A","BPQ020",
            "LBDTRSI","BPQ090D","LBDHDDSI",
            "MCQ160C","MCQ160D","MCQ160E","MCQ220",
            "DMDMARTZ","DMDEDUC2","SMQ020","SMQ040",
            "LBXSCR","LBDSCRSI","KIQ022","KIQ025",
            "LUXCAPM","LUXSMED",
            # ALNR inputs (EXACT set)
            "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI",
            "weekly_alcohol_grams"
        )
        data <- ensure_cols(data, need)
        
        data %>%
            mutate(
                # outcomes definition inputs
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                
                # cardiometabolic risks (NA-aware)
                risk_bmi_waist = if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                    NA,
                    (if_else(is.na(BMXBMI), FALSE, BMXBMI >= 25)) |
                        (RIAGENDR == 1 & if_else(is.na(BMXWAIST), FALSE, BMXWAIST >= 94)) |
                        (RIAGENDR == 2 & if_else(is.na(BMXWAIST), FALSE, BMXWAIST >= 80))
                ),
                
                risk_glucose_diabetes = if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7,9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7,9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7,9)),
                    NA,
                    (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
                        (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                
                any_bp_high = (
                    if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)
                ) | (
                    if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)
                ) | (
                    if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)
                ),
                
                risk_blood_pressure = if_else(
                    is.na(any_bp_high) & (is.na(BPQ040A) | BPQ040A %in% c(7,9)),
                    NA,
                    (if_else(is.na(any_bp_high), FALSE, any_bp_high)) |
                        (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                
                risk_triglycerides = if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9)),
                    NA,
                    (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                risk_low_hdl = if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9))),
                    NA,
                    ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                # ALNR (EXACT inputs + specified priority)
                Albumin_g_dl = get_albumin_gdl(cur_data()),
                NLR          = get_nlr(cur_data()),
                ALNR         = if_else(!is.na(Albumin_g_dl) & !is.na(NLR) & NLR > 0, Albumin_g_dl / NLR, NA_real_),
                ALNR_scaled  = as.vector(scale(ALNR))
            ) %>%
            rowwise() %>%
            mutate(num_risks = sum(c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl), na.rm = TRUE)) %>%
            ungroup() %>%
            mutate(
                has_one_plus_cardiometabolic_risk = if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA, num_risks > 0
                ),
                
                is_light_drinker_final = case_when(
                    is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
                    (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
                    TRUE ~ FALSE
                ),
                
                masld_group = case_when(
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
                    TRUE ~ NA_character_
                ),
                
                o1 = factor(if_else(masld_group == "MASLD", 1, 0, missing = 0),
                            levels = c(0,1), labels = c("No MASLD","MASLD")),
                
                slf_group = case_when(
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8 ~ "SLF",
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED < 8  ~ "non-SLF",
                    TRUE ~ NA_character_
                ),
                
                o2 = factor(if_else(slf_group == "SLF", 1, 0, missing = 0),
                            levels = c(0,1), labels = c("non-SLF","SLF")),
                
                # covariates
                age = RIDAGEYR,
                age_group = factor(if_else(RIDAGEYR >= 65, "≥65", "<65"), levels = c("<65","≥65")),
                gender = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
                
                race = fct_relevel(factor(case_when(
                    RIDRETH3 == 1 ~ "Mexican American",
                    RIDRETH3 == 2 ~ "Other Hispanic",
                    RIDRETH3 == 3 ~ "Non-Hispanic White",
                    RIDRETH3 == 4 ~ "Non-Hispanic Black",
                    RIDRETH3 == 6 ~ "Non-Hispanic Asian",
                    RIDRETH3 == 7 ~ "Other Race/Multi-Racial",
                    TRUE ~ NA_character_
                )), "Non-Hispanic White"),
                
                bmi_group = factor(if_else(BMXBMI < 25, "<25", "≥25"), levels = c("<25","≥25")),
                obesity   = factor(if_else(BMXBMI >= 30, "Yes", "No", missing = "No"), levels = c("No","Yes")),
                
                diabetes_mellitus = factor(
                    if_else(DIQ010==1 | DIQ050==1 | DIQ070==1 | (LBXGH>=6.5) | (LBXGLU>=126),
                            "Yes","No", missing="No"),
                    levels=c("No","Yes")
                ),
                
                hypertension = factor(
                    if_else(BPQ020==1 | BPQ040A==1, "Yes","No", missing="No"),
                    levels=c("No","Yes")
                ),
                
                cardiovascular_disease = factor(
                    if_else(MCQ160C==1 | MCQ160D==1 | MCQ160E==1, "Yes","No", missing="No"),
                    levels=c("No","Yes")
                ),
                
                cancer = factor(if_else(MCQ220==1, "Yes","No", missing="No"), levels=c("No","Yes")),
                
                marital_status = factor(case_when(
                    DMDMARTZ %in% c(1,6) ~ "Married/Living with partner",
                    DMDMARTZ %in% c(2,3,4) ~ "Widowed/Divorced/Separated",
                    DMDMARTZ == 5 ~ "Never married",
                    TRUE ~ NA_character_
                )),
                
                smoking_status = factor(case_when(
                    SMQ020 == 2 ~ "Non-smoker",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                    SMQ020 == 1 & SMQ040 %in% c(1,2) ~ "Current smoker",
                    TRUE ~ NA_character_
                ), levels=c("Non-smoker","Former smoker","Current smoker")),
                
                education = factor(case_when(
                    DMDEDUC2 == 1 ~ "Less than 9th grade",
                    DMDEDUC2 == 2 ~ "9-11th grade",
                    DMDEDUC2 == 3 ~ "High school graduate",
                    DMDEDUC2 == 4 ~ "Some college or AA degree",
                    DMDEDUC2 == 5 ~ "College graduate or above",
                    TRUE ~ NA_character_
                )),
                
                eGFR = {
                    creat <- case_when(
                        !is.na(LBXSCR) ~ LBXSCR,
                        !is.na(LBDSCRSI) ~ LBDSCRSI/88.4,
                        TRUE ~ NA_real_
                    )
                    if_else(!is.na(creat) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                            175*(creat^-1.154)*(RIDAGEYR^-0.203)*
                                if_else(RIAGENDR==2,0.742,1)*
                                if_else(RIDRETH1==4,1.212,1),
                            NA_real_)
                },
                
                chronic_kidney_disease = factor(
                    case_when(
                        !is.na(KIQ022) & KIQ022==1 ~ "Yes",
                        !is.na(KIQ025) & KIQ025==1 ~ "Yes",
                        !is.na(eGFR) & eGFR < 60 ~ "Yes",
                        TRUE ~ "No"
                    ),
                    levels=c("No","Yes")
                )
            ) %>%
            transmute(
                o1,o2,ALNR_scaled,
                age,age_group,gender,race,
                bmi_group,obesity,diabetes_mellitus,hypertension,
                cardiovascular_disease,chronic_kidney_disease,
                smoking_status,marital_status,education,cancer
            )
    }
    
    # ----------------------------- 3) MODELS / SUBGROUPS -----------------------------
    nice_lab <- c(
        age_group="Age", gender="Gender", bmi_group="BMI",
        hypertension="Hypertension", diabetes_mellitus="Diabetes",
        cardiovascular_disease="CVD", chronic_kidney_disease="Renal disease",
        obesity="Obesity", smoking_status="Smoking", race="Race",
        marital_status="Marital status", education="Education", cancer="Cancer"
    )
    
    covariates_model <- c(
        "age","gender","race","marital_status","education","smoking_status",
        "bmi_group","obesity","diabetes_mellitus","hypertension","cancer",
        "chronic_kidney_disease","cardiovascular_disease"
    )
    
    subgroup_vars <- c(
        "age_group","gender","bmi_group","hypertension","diabetes_mellitus",
        "cardiovascular_disease","chronic_kidney_disease","obesity",
        "smoking_status","race","marital_status","education","cancer"
    )
    
    safe_fit <- function(formula, df) {
        ok <- TRUE
        m <- try(suppressWarnings(glm(formula, data=df, family=binomial())), silent=TRUE)
        if (inherits(m,"try-error") ||
            any(!is.finite(coef(m))) ||
            any(!is.finite(diag(vcov(m))))) ok <- FALSE
        
        if (ok) {
            est <- unname(coef(m)["ALNR_scaled"])
            se  <- sqrt(vcov(m)["ALNR_scaled","ALNR_scaled"])
            return(list(
                or = exp(est),
                lo = exp(est - 1.96*se),
                hi = exp(est + 1.96*se),
                p  = summary(m)$coefficients["ALNR_scaled","Pr(>|z|)"]
            ))
        }
        
        mf <- suppressWarnings(try(logistf(formula, data=df), silent=TRUE))
        if (inherits(mf,"try-error")) return(list(or=NA, lo=NA, hi=NA, p=NA))
        
        list(
            or = exp(mf$coefficients["ALNR_scaled"]),
            lo = exp(mf$ci.lower["ALNR_scaled"]),
            hi = exp(mf$ci.upper["ALNR_scaled"]),
            p  = mf$prob["ALNR_scaled"]
        )
    }
    
    interaction_p <- function(df, var, ycol="y") {
        if (!var %in% names(df) || n_distinct(df[[var]]) < 2) return(NA_real_)
        covs <- setdiff(covariates_model, var)
        covs <- covs[sapply(df[, covs, drop=FALSE], n_distinct) >= 2]
        
        f0 <- if (length(covs)==0) as.formula(paste(ycol, "~ ALNR_scaled +", var))
        else as.formula(paste(ycol, "~ ALNR_scaled +", paste(covs, collapse=" + "), "+", var))
        
        f1 <- if (length(covs)==0) as.formula(paste(ycol, "~ ALNR_scaled *", var))
        else as.formula(paste(ycol, "~ ALNR_scaled *", var, "+", paste(covs, collapse=" + ")))
        
        m0 <- try(suppressWarnings(glm(f0, data=df, family=binomial())), silent=TRUE)
        m1 <- try(suppressWarnings(glm(f1, data=df, family=binomial())), silent=TRUE)
        if (inherits(m0,"try-error") || inherits(m1,"try-error")) return(NA_real_)
        
        ll0 <- as.numeric(logLik(m0)); ll1 <- as.numeric(logLik(m1))
        dfdeg <- length(coef(m1)) - length(coef(m0))
        pchisq(2*(ll1-ll0), df=dfdeg, lower.tail=FALSE)
    }
    
    prepare_result <- function(df, outcome_factor, positive_label) {
        
        df <- df %>%
            drop_na(ALNR_scaled) %>%
            filter(!is.na(.data[[outcome_factor]])) %>%
            mutate(y = as.integer(.data[[outcome_factor]] == positive_label))
        
        covs_overall <- covariates_model[
            sapply(df[, covariates_model, drop=FALSE], n_distinct) >= 2
        ]
        
        f_overall <- if (length(covs_overall)==0) {
            as.formula("y ~ ALNR_scaled")
        } else {
            as.formula(paste("y ~ ALNR_scaled +", paste(covs_overall, collapse=" + ")))
        }
        
        ov <- safe_fit(f_overall, df)
        
        res_overall <- tibble(
            group="Overall",
            group_label="Overall",
            level="Overall",
            n=nrow(df),
            or=ov$or, lo=ov$lo, hi=ov$hi, p=ov$p,
            p_int=NA_real_
        )
        
        sgv <- intersect(subgroup_vars, names(df))
        
        sub <- map_dfr(sgv, function(var) {
            lvls <- levels(df[[var]])
            pint <- interaction_p(df, var, "y")
            
            rows <- map_dfr(lvls, function(lv) {
                dd <- df %>% filter(.data[[var]] == lv)
                
                if (nrow(dd) < 20 || n_distinct(dd$y) < 2)
                    return(tibble(or=NA, lo=NA, hi=NA, p=NA, n=nrow(dd)))
                
                covs <- setdiff(covs_overall, var)
                covs <- covs[sapply(dd[, covs, drop=FALSE], n_distinct) >= 2]
                
                f <- if (length(covs)==0) as.formula("y ~ ALNR_scaled")
                else as.formula(paste("y ~ ALNR_scaled +", paste(covs, collapse=" + ")))
                
                ft <- safe_fit(f, dd)
                tibble(or=ft$or, lo=ft$lo, hi=ft$hi, p=ft$p, n=nrow(dd))
            })
            
            rows %>%
                mutate(
                    group = var,
                    group_label = unname(nice_lab[var]),
                    level = lvls,
                    p_int = c(pint, rep(NA_real_, length(lvls)-1))
                )
        })
        
        bind_rows(res_overall, sub)
    }
    
    # ----------------------------- 4) 4-PANEL PLOT (NO N COLUMN) + AUTO SAVE -----------------------------
    panel_metrics <- function(res){
        max_level <- max(nchar(as.character(res$level)), na.rm=TRUE)
        max_level <- ifelse(is.finite(max_level), max_level, 30)
        
        w_name <- 1.70 + 0.040 * max(0, max_level - 28)
        w_mid  <- 2.10
        w_aor  <- 1.15
        w_p    <- 0.90
        
        width_in  <- 10.2 + 0.10 * max(0, max_level - 28)
        width_in  <- max(10.2, min(15.0, width_in))
        height_in <- max(5.0, 0.235 * nrow(res) + 1.1)
        
        list(widths=c(w_name,w_mid,w_aor,w_p), w=width_in, h=height_in)
    }
    
    plot_forest_4panel <- function(res){
        
        txt_body   <- 5.1
        txt_header <- 5.3
        base_mid   <- 12
        
        res <- res %>%
            mutate(group_label=factor(group_label, levels=c("Overall", setdiff(unique(group_label),"Overall")))) %>%
            group_by(group_label) %>% mutate(row_in_group=row_number()) %>%
            ungroup() %>%
            mutate(row_id = row_number(),
                   y = rev(seq_len(n())))
        
        aor_txt <- function(or,lo,hi){
            if (any(is.na(c(or,lo,hi)))) return("")
            paste0(formatC(or,digits=2,format="f")," (",
                   formatC(lo,digits=2,format="f"),"–",
                   formatC(hi,digits=2,format="f"),")")
        }
        res$aor_label <- pmap_chr(list(res$or,res$lo,res$hi), ~aor_txt(..1,..2,..3))
        y_top <- max(res$y) + 1.25
        
        # row shading (for all panels)
        shade <- res %>%
            mutate(
                ymin = y - 0.5,
                ymax = y + 0.5,
                fill_row = if_else(row_id %% 2L == 0L, "grey95", "white")
            )
        
        # group separators (skip first header row)
        sep_y <- res %>% filter(row_in_group==1) %>% pull(y)
        sep_y <- sep_y[sep_y != max(res$y)] + 0.5
        
        # Panel 1: subgroup names
        x_group_end   <- 1.00
        x_level_start <- 1.10
        x_name_xlim   <- c(0, 3.8)
        
        p_name <- ggplot() +
            geom_rect(
                data=shade,
                aes(xmin=x_name_xlim[1], xmax=x_name_xlim[2], ymin=ymin, ymax=ymax),
                fill=shade$fill_row, color=NA
            ) +
            geom_hline(yintercept = sep_y, linewidth=0.25, color="grey80") +
            geom_text(
                data=res,
                aes(y=y, x=x_group_end, label=if_else(row_in_group==1, as.character(group_label), "")),
                hjust=1, fontface="bold", size=txt_body, lineheight=0.95
            ) +
            geom_text(
                data=res,
                aes(y=y, x=x_level_start, label=if_else(group_label=="Overall", "", as.character(level))),
                hjust=0, size=txt_body, lineheight=0.95
            ) +
            annotate("text", x=x_name_xlim[1], y=y_top, label="Subgroups", hjust=0, fontface="bold", size=txt_header) +
            scale_x_continuous(breaks=NULL) +
            scale_y_continuous(breaks=NULL) +
            coord_cartesian(xlim=x_name_xlim, clip="off") +
            theme_void() +
            theme(plot.margin=margin(2,6,2,8))
        
        # Panel 2: forest
        xlo <- max(0.2, min(res$lo, na.rm=TRUE) * 0.80)
        xhi <- min(8.0, max(res$hi, na.rm=TRUE) * 1.20)
        
        p_mid <- ggplot() +
            geom_rect(
                data=shade,
                aes(xmin=xlo, xmax=xhi, ymin=ymin, ymax=ymax),
                fill=shade$fill_row, color=NA
            ) +
            geom_hline(yintercept = sep_y, linewidth=0.25, color="grey80") +
            geom_vline(xintercept=1, linetype="dashed", linewidth=0.5, color="grey30") +
            geom_errorbarh(
                data=res,
                aes(xmin=lo, xmax=hi, y=y),
                height=0.24, linewidth=0.65, na.rm=TRUE, color="grey20"
            ) +
            geom_point(
                data=res,
                aes(x=or, y=y),
                size=2.25, na.rm=TRUE, color="black"
            ) +
            scale_x_log10(
                limits=c(xlo, xhi),
                breaks=c(0.5,0.7,1,1.5,2,3,5),
                labels=c("0.5","0.7","1","1.5","2","3","5"),
                expand=expansion(mult=c(0,0))
            ) +
            scale_y_continuous(breaks=NULL, expand=expansion(mult=c(0.01,0.01))) +
            labs(x="Adjusted OR (95% CI) per 1 SD increase in ALNR", y=NULL) +
            theme_minimal(base_size=base_mid) +
            theme(
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_text(face="bold"),
                axis.text.x = element_text(color="black"),
                plot.margin = margin(2,6,2,6)
            )
        
        # Panel 3: aOR text
        x_aor_xlim <- c(0, 1)
        
        p_aor <- ggplot() +
            geom_rect(
                data=shade,
                aes(xmin=x_aor_xlim[1], xmax=x_aor_xlim[2], ymin=ymin, ymax=ymax),
                fill=shade$fill_row, color=NA
            ) +
            geom_hline(yintercept = sep_y, linewidth=0.25, color="grey80") +
            geom_text(
                data=res,
                aes(y=y, x=0, label=aor_label),
                hjust=0, size=txt_body, lineheight=0.95, color="black"
            ) +
            annotate("text", x=0, y=y_top, label="aOR (95% CI)", hjust=0, fontface="bold", size=txt_header) +
            scale_x_continuous(breaks=NULL) +
            scale_y_continuous(breaks=NULL) +
            coord_cartesian(xlim=x_aor_xlim, clip="off") +
            theme_void() +
            theme(plot.margin=margin(2,6,2,4))
        
        # Panel 4: P interaction
        x_p_xlim <- c(0, 1)
        
        p_p <- ggplot() +
            geom_rect(
                data=shade,
                aes(xmin=x_p_xlim[1], xmax=x_p_xlim[2], ymin=ymin, ymax=ymax),
                fill=shade$fill_row, color=NA
            ) +
            geom_hline(yintercept = sep_y, linewidth=0.25, color="grey80") +
            geom_text(
                data=res,
                aes(y=y, x=0, label=ifelse(row_in_group==1 & !is.na(p_int), scales::pvalue(p_int, accuracy=0.001), "")),
                hjust=0, size=txt_body, lineheight=0.95, color="black"
            ) +
            annotate("text", x=0, y=y_top, label="P for interaction", hjust=0, fontface="bold", size=txt_header) +
            scale_x_continuous(breaks=NULL) +
            scale_y_continuous(breaks=NULL) +
            coord_cartesian(xlim=x_p_xlim, clip="off") +
            theme_void() +
            theme(plot.margin=margin(2,2,2,2))
        
        ws <- panel_metrics(res)$widths
        (p_name + p_mid + p_aor + p_p) + patchwork::plot_layout(widths = ws)
    }
    
    save_forest_4panel <- function(p, res, file, dpi=300){
        m <- panel_metrics(res)
        ggsave(file, p, width=m$w, height=m$h, units="in", dpi=dpi, bg="white", limitsize=FALSE)
    }
    
    # ----------------------------- 5) RUN (NO SPLIT) -----------------------------
    if (!exists("data")) stop("Thiếu object `data` trong Environment.")
    
    dir.create("plots", showWarnings = FALSE)
    
    data_cleaned  <- apply_exclusion_criteria(data)
    data_featured <- create_features(data_cleaned)
    
    # MASLD (full cohort)
    res_masld <- prepare_result(data_featured, "o1", "MASLD")
    p_masld   <- plot_forest_4panel(res_masld)
    print(p_masld)
    save_forest_4panel(p_masld, res_masld, file.path("plots","forest_MASLD_ALNR.png"), dpi=300)
    
    # SLF (MASLD only)
    res_slf <- prepare_result(data_featured %>% filter(o1=="MASLD"), "o2", "SLF")
    p_slf   <- plot_forest_4panel(res_slf)
    print(p_slf)
    save_forest_4panel(p_slf, res_slf, file.path("plots","forest_SLF_ALNR.png"), dpi=300)
    
    ```
    
- **Medication**
    
    ```r
    # ==============================================================================
    # END-TO-END (UPDATED AS REQUESTED)
    # 1) NO train/validation split (use FULL cleaned dataset)
    # 2) Replace ALL CURI/hs-CRP/uric acid logic with Albumin/NLR
    # 3) Exclusion "key measurements" use ONLY Alb/NLR inputs (plus FibroScan for outcomes):
    #    - LBXLYPCT  % lymphocyte
    #    - LBXNEPCT  % neutrophil
    #    - LBDLYMNO  lymphocyte (absolute)
    #    - LBDNENO   neutrophil (absolute)
    #    - LBXSAL    Albumin (g/dL)
    #    - LBDSALSI  Albumin (g/L)
    #    NLR priority: abs (LBDNENO/LBDLYMNO), fallback % (LBXNEPCT/LBXLYPCT)
    #    Albumin priority: LBXSAL, fallback LBDSALSI/10
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(tidyverse)
        library(mediation)
    })
    
    # ----------------------------- 0) HELPERS -----------------------------
    ensure_cols <- function(df, cols){
        miss <- setdiff(cols, names(df))
        if (length(miss)) for (m in miss) df[[m]] <- NA
        df
    }
    num_safe <- function(x) suppressWarnings(as.numeric(x))
    
    get_albumin_gdl_vec <- function(df){
        lbx <- num_safe(df$LBXSAL)      # g/dL
        lbs <- num_safe(df$LBDSALSI)    # g/L
        ifelse(!is.na(lbx), lbx,
               ifelse(!is.na(lbs), lbs/10, NA_real_))
    }
    
    get_nlr_vec <- function(df){
        neu_abs <- num_safe(df$LBDNENO)
        lym_abs <- num_safe(df$LBDLYMNO)
        neu_pct <- num_safe(df$LBXNEPCT)
        lym_pct <- num_safe(df$LBXLYPCT)
        
        nlr_abs <- ifelse(!is.na(neu_abs) & !is.na(lym_abs) & lym_abs > 0, neu_abs/lym_abs, NA_real_)
        nlr_pct <- ifelse(!is.na(neu_pct) & !is.na(lym_pct) & lym_pct > 0, neu_pct/lym_pct, NA_real_)
        
        ifelse(!is.na(nlr_abs), nlr_abs, nlr_pct)
    }
    
    # ----------------------------- PART 1: EXCLUSION -----------------------------
    cat("--- PART 1: Starting initial data filtering process ---\n")
    initial_rows_total <- nrow(data)
    
    # Remove completely empty columns
    data_no_empty_cols <- data[, colSums(!is.na(data)) > 0, drop = FALSE]
    
    # Ensure columns that will be referenced in filtering
    data_no_empty_cols <- ensure_cols(data_no_empty_cols, c(
        # hepatitis (optional)
        "RHD143","LBXHBC","LBDHBG","LBXHCR",
        # FibroScan for outcomes
        "LUXSMED","LUXCAPM",
        # Alb/NLR inputs (EXACT set)
        "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI",
        # alcohol + demographics (optional but used)
        "ALQ121","ALQ130","RIAGENDR","RIDAGEYR"
    ))
    
    # 1.1 Exclude viral hepatitis (only if any hep var exists with non-NA data)
    hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR")
    hep_vars_exist <- intersect(hep_vars, names(data_no_empty_cols))
    rows_before <- nrow(data_no_empty_cols)
    
    if (length(hep_vars_exist) > 0) {
        data_after_hep_filter <- data_no_empty_cols %>%
            filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
    } else {
        data_after_hep_filter <- data_no_empty_cols
    }
    
    cat("-> Removed", rows_before - nrow(data_after_hep_filter),
        "subjects with evidence of viral hepatitis.\n")
    
    # 1.2 Exclude missing core liver measurement variables (needed for outcomes)
    vars_must_exist <- c("LUXSMED", "LUXCAPM")
    missing_vars <- setdiff(vars_must_exist, names(data_after_hep_filter))
    if (length(missing_vars) > 0) {
        stop(paste("Error: Required columns for analysis not found:", paste(missing_vars, collapse = ", ")))
    }
    
    rows_before <- nrow(data_after_hep_filter)
    data_after_na_drop <- data_after_hep_filter %>%
        drop_na(all_of(vars_must_exist))
    
    cat("-> Removed", rows_before - nrow(data_after_na_drop),
        "rows with missing core measurement data (LSM, CAP).\n")
    
    # 1.3 Exclude heavy alcohol drinkers (skip if ALQ missing)
    data_with_alcohol <- data_after_na_drop %>%
        mutate(
            days_per_week_alcohol = case_when(
                is.na(ALQ121) ~ NA_real_,
                ALQ121 == 0 ~ 0,
                ALQ121 == 1 ~ 7,
                ALQ121 == 2 ~ 6,
                ALQ121 == 3 ~ 3.5,
                ALQ121 == 4 ~ 2,
                ALQ121 == 5 ~ 1,
                ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                ALQ121 == 7 ~ 1 / (30.4375 / 7),
                ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                TRUE ~ NA_real_
            ),
            ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, as.numeric(ALQ130)),
            weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
        )
    
    rows_before <- nrow(data_with_alcohol)
    data_after_alcohol_filter <- data_with_alcohol %>%
        filter(
            is.na(weekly_alcohol_grams) |
                (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
        )
    
    cat("-> Removed", rows_before - nrow(data_after_alcohol_filter),
        "rows for exceeding alcohol consumption thresholds.\n")
    
    # 1.4 Exclude subjects under 18
    rows_before <- nrow(data_after_alcohol_filter)
    data_age_ok <- data_after_alcohol_filter %>%
        filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
    
    cat("-> Removed", rows_before - nrow(data_age_ok), "rows for being under 18 years old.\n")
    
    # 1.5 Key measurements exclusion: ONLY Alb/NLR inputs (compute ALB_NLR and require non-missing)
    #     Albumin priority: LBXSAL else LBDSALSI/10
    #     NLR priority: abs else %
    rows_before <- nrow(data_age_ok)
    
    data_cleaned <- data_age_ok %>%
        mutate(
            albumin_g_dl = get_albumin_gdl_vec(cur_data()),
            NLR = get_nlr_vec(cur_data()),
            ALB_NLR = if_else(!is.na(albumin_g_dl) & !is.na(NLR) & NLR > 0, albumin_g_dl / NLR, NA_real_)
        ) %>%
        filter(!is.na(ALB_NLR))
    
    cat("-> Removed", rows_before - nrow(data_cleaned),
        "rows with insufficient Alb/NLR inputs (cannot compute Albumin/NLR).\n")
    
    cat("--- Data filtering process complete. Total rows removed:",
        initial_rows_total - nrow(data_cleaned),
        ". Remaining rows:", nrow(data_cleaned), "---\n\n")
    
    # ----------------------------- PART 2: FEATURE ENGINEERING -----------------------------
    cat("--- PART 2: Starting robust feature engineering for analysis ---\n")
    
    # Ensure ALL referenced variables exist to avoid hard errors (non-key vars can be NA)
    data_cleaned <- ensure_cols(data_cleaned, c(
        # glucose / lipids / insulin (for mediators)
        "LBXGLU","LBDGLUSI","LBXTR","LBDTRSI","LBDHDD","LBDHDDSI","LBXIN",
        # cardiometabolic / anthropometrics
        "BMXBMI","BMXWAIST","LBXGH","DIQ010","DIQ050","DIQ070",
        "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3","BPQ040A","BPQ020","BPQ090D",
        # covariates
        "LBXSCR","LBDSCRSI","KIQ022","KIQ025",
        "MCQ160C","MCQ160D","MCQ160E","MCQ220",
        "RIDRETH1","RIDRETH3","DMDMARTZ","DMDEDUC2","SMQ020","SMQ040"
    ))
    
    data_featured <- data_cleaned %>%
        mutate(
            # A) Biomarkers for mediators (NA-safe conversions)
            glucose_mg_dl = case_when(
                !is.na(LBXGLU) ~ as.numeric(LBXGLU),
                !is.na(LBDGLUSI) ~ as.numeric(LBDGLUSI) / 0.05551,
                TRUE ~ NA_real_
            ),
            triglycerides_mg_dl = case_when(
                !is.na(LBXTR) ~ as.numeric(LBXTR),
                !is.na(LBDTRSI) ~ as.numeric(LBDTRSI) * 88.57,
                TRUE ~ NA_real_
            ),
            hdl_mg_dl = case_when(
                !is.na(LBDHDD) ~ as.numeric(LBDHDD),
                !is.na(LBDHDDSI) ~ as.numeric(LBDHDDSI) * 38.67,
                TRUE ~ NA_real_
            ),
            
            # B) Predictor = Albumin/NLR (already computed; recompute deterministically to keep logic local)
            albumin_g_dl = get_albumin_gdl_vec(cur_data()),
            NLR          = get_nlr_vec(cur_data()),
            ALB_NLR      = if_else(!is.na(albumin_g_dl) & !is.na(NLR) & NLR > 0, albumin_g_dl / NLR, NA_real_),
            
            # C) Mediators (robust tryCatch)
            METS_IR = tryCatch({
                (log(2 * glucose_mg_dl + triglycerides_mg_dl) * BMXBMI) / log(hdl_mg_dl)
            }, error = function(e) NA_real_),
            HOMA_IR = tryCatch({
                (glucose_mg_dl * LBXIN) / 405
            }, error = function(e) NA_real_),
            
            # D) Outcomes (MASLD, SLF)
            has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
            
            has_one_plus_cardiometabolic_risk = {
                risk_bmi_waist <- (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)
                risk_glucose_diabetes <- (glucose_mg_dl >= 100) | (LBXGH >= 5.7) | (DIQ010 == 1) | (DIQ050 == 1) | (DIQ070 == 1)
                any_bp_high <- (BPXOSY1 >= 130 | BPXODI1 >= 85) | (BPXOSY2 >= 130 | BPXODI2 >= 85) | (BPXOSY3 >= 130 | BPXODI3 >= 85)
                risk_blood_pressure <- any_bp_high | (BPQ040A == 1)
                risk_triglycerides <- (triglycerides_mg_dl >= 150) | (BPQ090D == 1)
                risk_low_hdl <- ((RIAGENDR == 1 & hdl_mg_dl < 40) | (RIAGENDR == 2 & hdl_mg_dl < 50)) | (BPQ090D == 1)
                num_risks <- rowSums(cbind(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
                                     na.rm = TRUE)
                num_risks > 0
            },
            
            masld_group = case_when(
                has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
                has_hepatic_steatosis == FALSE ~ "non-MASLD",
                has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                TRUE ~ NA_character_
            ),
            masld_numeric = if_else(masld_group == "MASLD", 1, 0),
            
            slf_group = case_when(
                masld_group == "MASLD" & LUXSMED >= 8.0 ~ "SLF",
                masld_group == "MASLD" & LUXSMED < 8.0 ~ "non-SLF",
                TRUE ~ NA_character_
            ),
            slf_binary = if_else(slf_group == "SLF", 1, 0),
            
            # E) Covariates
            eGFR = {
                creatinine_mg_dl <- case_when(
                    !is.na(LBXSCR) ~ as.numeric(LBXSCR),
                    !is.na(LBDSCRSI) ~ as.numeric(LBDSCRSI) / 88.4,
                    TRUE ~ NA_real_
                )
                if_else(!is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                        175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                            if_else(RIAGENDR == 2, 0.742, 1) *
                            if_else(RIDRETH1 == 4, 1.212, 1),
                        NA_real_)
            },
            age = RIDAGEYR,
            gender = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
            race = factor(case_when(
                RIDRETH1 == 1 ~ "Mexican American",
                RIDRETH1 == 2 ~ "Other Hispanic",
                RIDRETH1 == 3 ~ "Non-Hispanic White",
                RIDRETH1 == 4 ~ "Non-Hispanic Black",
                RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
                TRUE ~ "Other Race/Multi-Racial"
            )),
            marital_status = factor(case_when(
                DMDMARTZ %in% c(1, 6) ~ "Married/Living with partner",
                DMDMARTZ %in% c(2, 3, 4) ~ "Widowed/Divorced/Separated",
                DMDMARTZ == 5 ~ "Never married",
                TRUE ~ NA_character_
            )),
            education = factor(case_when(
                DMDEDUC2 == 1 ~ "Less than 9th grade",
                DMDEDUC2 == 2 ~ "9-11th grade",
                DMDEDUC2 == 3 ~ "High school graduate",
                DMDEDUC2 == 4 ~ "Some college or AA degree",
                DMDEDUC2 == 5 ~ "College graduate or above",
                TRUE ~ NA_character_
            )),
            smoking_status = factor(case_when(
                SMQ020 == 2 ~ "Non-smoker",
                SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
                TRUE ~ NA_character_
            ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
            obesity = factor(if_else(BMXBMI >= 30, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
            diabetes_mellitus = factor(if_else(
                DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                    (!is.na(LBXGH) & LBXGH >= 6.5) |
                    (!is.na(glucose_mg_dl) & glucose_mg_dl >= 126),
                "Yes", "No", missing = "No"
            ), levels = c("No", "Yes")),
            hypertension = factor(if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"),
                                  levels = c("No", "Yes")),
            chronic_kidney_disease = factor(if_else(
                (!is.na(KIQ022) & KIQ022 == 1) | (!is.na(KIQ025) & KIQ025 == 1) | (!is.na(eGFR) & eGFR < 60),
                "Yes", "No", missing = "No"
            ), levels = c("No", "Yes")),
            cancer = factor(if_else(MCQ220 == 1, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
            cardiovascular_disease = factor(if_else(MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1, "Yes", "No", missing = "No"),
                                            levels = c("No", "Yes"))
        )
    
    cat("--- Feature engineering complete ---\n\n")
    
    # Enforce predictor availability (should already hold from exclusion, but keep as final guard)
    rows_before <- nrow(data_featured)
    data_for_full_analysis <- data_featured %>% filter(!is.na(ALB_NLR))
    cat("-> Removed", rows_before - nrow(data_for_full_analysis),
        "rows with missing Albumin/NLR (predictor).\n")
    cat("-> Rows available after cleaning + predictor availability:", nrow(data_for_full_analysis), "\n\n")
    
    # ----------------------------- PART 3: MEDIATION ANALYSIS -----------------------------
    cat("--- PART 3: Performing Mediation Analysis on the FULL dataset ---\n")
    
    outcomes_for_analysis  <- c("masld_numeric", "slf_binary")
    mediators_for_analysis <- c("METS_IR", "HOMA_IR")
    mediation_results <- list()
    
    covariates_list_full <- c(
        "age", "gender", "race", "marital_status", "education", "smoking_status",
        "obesity", "diabetes_mellitus", "hypertension", "cancer",
        "chronic_kidney_disease", "cardiovascular_disease"
    )
    
    for (outcome_var in outcomes_for_analysis) {
        mediation_results[[outcome_var]] <- list()
        
        for (mediator_var in mediators_for_analysis) {
            
            cat(sprintf("\n------------------------------------------------------------\n"))
            cat(sprintf(">>> Starting analysis for Outcome: %s & Mediator: %s <<<\n",
                        toupper(sub("_.*", "", outcome_var)), mediator_var))
            
            vars_to_exclude <- character()
            if (mediator_var == "METS_IR") {
                vars_to_exclude <- c("obesity", "diabetes_mellitus", "hypertension")
            } else if (mediator_var == "HOMA_IR") {
                vars_to_exclude <- "diabetes_mellitus"
            }
            covariates_to_use <- setdiff(covariates_list_full, vars_to_exclude)
            
            if (length(vars_to_exclude) > 0) {
                cat(sprintf("       NOTE: Excluding '%s' from covariates for the %s model to prevent multicollinearity.\n",
                            paste(vars_to_exclude, collapse = ", "), mediator_var))
            }
            
            full_vars <- c("ALB_NLR", outcome_var, mediator_var, covariates_to_use)
            data_subset <- data_for_full_analysis %>%
                dplyr::select(any_of(full_vars)) %>%
                na.omit()
            
            actual_n <- nrow(data_subset)
            cat(sprintf("       Number of complete observations for this analysis: %d\n", actual_n))
            
            if (actual_n < 50) {
                cat(sprintf("--- WARNING: Insufficient data (n=%d) to run analysis. Skipping.\n", actual_n))
                next
            }
            
            valid_covariates <- character()
            for (cov in covariates_to_use) {
                if (cov %in% names(data_subset)) {
                    if (!is.factor(data_subset[[cov]]) || dplyr::n_distinct(data_subset[[cov]]) >= 2) {
                        valid_covariates <- c(valid_covariates, cov)
                    }
                }
            }
            cov_part <- paste(valid_covariates, collapse = " + ")
            
            cat("       Running regression models and the mediate() algorithm...\n")
            
            model.m <- NULL
            model.y <- NULL
            
            analysis_result <- tryCatch({
                model_m_formula <- as.formula(paste(mediator_var, "~ ALB_NLR +", cov_part))
                model.m <- lm(model_m_formula, data = data_subset)
                
                model_y_formula <- as.formula(paste(outcome_var, "~ ALB_NLR +", mediator_var, "+", cov_part))
                model.y <- glm(model_y_formula, data = data_subset, family = binomial(link = "logit"))
                
                mediate(model.m, model.y, treat = "ALB_NLR", mediator = mediator_var, boot = TRUE, sims = 500)
            }, error = function(e) {
                cat(sprintf("\n!!! ERROR running mediate() for '%s' with outcome '%s': %s\n",
                            mediator_var, outcome_var, e$message))
                return(NULL)
            })
            
            if (!is.null(analysis_result)) {
                mediation_results[[outcome_var]][[mediator_var]] <- list(
                    result  = analysis_result,
                    model_m = model.m,
                    model_y = model.y
                )
                cat("       Analysis complete.\n")
            }
        }
    }
    
    cat("\n--- Completed all mediation analyses ---\n\n")
    
    # ----------------------------- PART 4: PLOT MEDIATION DIAGRAMS -----------------------------
    cat("--- PART 4: Plotting publication-quality mediation diagrams ---\n")
    
    plot_publication_mediation <- function(mediation_bundle,
                                           predictor_var, predictor_label,
                                           mediator_var, mediator_label,
                                           outcome_label, title_letter,
                                           color_palette) {
        
        s   <- summary(mediation_bundle$result)
        s_m <- summary(mediation_bundle$model_m)
        s_y <- summary(mediation_bundle$model_y)
        
        format_p <- function(p_val) {
            if (is.na(p_val)) return("p = NA")
            if (p_val < 0.001) return("p < 0.001")
            sprintf("p = %.3f", p_val)
        }
        
        # Path a (β)
        path_a_text <- sprintf("β = %.2f\n(%s)",
                               s_m$coefficients[predictor_var, "Estimate"],
                               format_p(s_m$coefficients[predictor_var, "Pr(>|t|)"]))
        
        # Path b (OR)
        coef_b   <- s_y$coefficients[mediator_var, "Estimate"]
        se_b     <- s_y$coefficients[mediator_var, "Std. Error"]
        or_b     <- exp(coef_b)
        ci_b_low <- exp(coef_b - 1.96 * se_b)
        ci_b_high<- exp(coef_b + 1.96 * se_b)
        path_b_text <- sprintf("OR = %.2f (95%% CI: %.2f-%.2f)\n%s",
                               or_b, ci_b_low, ci_b_high,
                               format_p(s_y$coefficients[mediator_var, "Pr(>|z|)"]))
        
        # Direct (ADE) (OR)
        coef_c   <- s_y$coefficients[predictor_var, "Estimate"]
        se_c     <- s_y$coefficients[predictor_var, "Std. Error"]
        or_c     <- exp(coef_c)
        ci_c_low <- exp(coef_c - 1.96 * se_c)
        ci_c_high<- exp(coef_c + 1.96 * se_c)
        ade_text <- sprintf("Direct Effect (ADE)\nOR = %.2f (95%% CI: %.2f-%.2f)\n%s",
                            or_c, ci_c_low, ci_c_high,
                            format_p(s_y$coefficients[predictor_var, "Pr(>|z|)"]))
        
        # Indirect (ACME) + Proportion mediated
        acme_text <- sprintf("Indirect Effect (ACME) = %.3f, %s\n95%% CI [%.3f, %.3f]",
                             s$z.avg, format_p(s$z.avg.p), s$z.avg.ci[1], s$z.avg.ci[2])
        prop_text <- sprintf("Proportion Mediated = %.1f%%\n95%% CI [%.1f%%, %.1f%%]",
                             s$n.avg * 100, s$n.avg.ci[1] * 100, s$n.avg.ci[2] * 100)
        
        par(family = "sans", bg = "white")
        plot(NULL, xlim = c(0, 20), ylim = c(0, 20), axes = FALSE, xlab = "", ylab = "", asp = 1)
        
        text(1, 19.5, title_letter, cex = 2.5, font = 2, adj = c(0, 1))
        
        box_h <- 1.8; box_w <- 4.5
        pos <- list(
            x = c(x = 5,  y = 5),   # predictor
            m = c(x = 10, y = 17),  # mediator
            y = c(x = 15, y = 5)    # outcome
        )
        
        rect(pos$x['x'] - box_w, pos$x['y'] - box_h, pos$x['x'] + box_w, pos$x['y'] + box_h,
             col = color_palette$predictor, border = "black", lwd = 2)
        text(pos$x['x'], pos$x['y'], predictor_label, cex = 1.5, font = 2)
        
        rect(pos$m['x'] - box_w, pos$m['y'] - box_h, pos$m['x'] + box_w, pos$m['y'] + box_h,
             col = color_palette$mediator, border = "black", lwd = 2)
        text(pos$m['x'], pos$m['y'], mediator_label, cex = 1.5, font = 2)
        
        rect(pos$y['x'] - box_w, pos$y['y'] - box_h, pos$y['x'] + box_w, pos$y['y'] + box_h,
             col = color_palette$outcome, border = "black", lwd = 2)
        text(pos$y['x'], pos$y['y'], outcome_label, cex = 1.5, font = 2)
        
        arrow_coords <- list(
            a = list(start = c(pos$x['x'], pos$x['y'] + box_h), end = c(pos$m['x'], pos$m['y'] - box_h)),
            b = list(start = c(pos$m['x'], pos$m['y'] - box_h), end = c(pos$y['x'], pos$y['y'] + box_h)),
            c = list(start = c(pos$x['x'] + box_w, pos$x['y']), end = c(pos$y['x'] - box_w, pos$y['y']))
        )
        
        place_text_on_arrow <- function(start, end, label, offset_mult = 1.8, cex_val = 1.1) {
            mid_point <- (start + end) / 2
            angle <- atan2(end[2] - start[2], end[1] - start[1]) * 180 / pi
            perp_vec <- c(-(end[2] - start[2]), end[1] - start[1])
            norm_perp <- perp_vec / sqrt(sum(perp_vec^2))
            text_pos <- mid_point + norm_perp * offset_mult
            text(text_pos[1], text_pos[2], label, srt = angle, cex = cex_val, font = 1,
                 adj = c(0.5, 0.5), xpd = TRUE)
        }
        
        arrow_lwd <- 2
        arrows(arrow_coords$a$start[1], arrow_coords$a$start[2], arrow_coords$a$end[1], arrow_coords$a$end[2],
               length = 0.1, lwd = arrow_lwd, angle = 25)
        place_text_on_arrow(arrow_coords$a$start, arrow_coords$a$end, path_a_text)
        
        arrows(arrow_coords$b$start[1], arrow_coords$b$start[2], arrow_coords$b$end[1], arrow_coords$b$end[2],
               length = 0.1, lwd = arrow_lwd, angle = 25)
        place_text_on_arrow(arrow_coords$b$start, arrow_coords$b$end, path_b_text)
        
        arrows(arrow_coords$c$start[1], arrow_coords$c$start[2], arrow_coords$c$end[1], arrow_coords$c$end[2],
               length = 0.1, lwd = arrow_lwd, angle = 25)
        text(mean(c(arrow_coords$c$start[1], arrow_coords$c$end[1])),
             arrow_coords$c$start[2] - 2.0, ade_text, cex = 1.1, font = 1, adj = c(0.5, 0.5))
        
        text(10, 11.5, acme_text, cex = 1.1, font = 2, adj = c(0.5, 0.5))
        text(10, 9.5,  prop_text, cex = 1.1, font = 1, adj = c(0.5, 0.5))
    }
    
    has_any <- any(sapply(mediation_results, function(x) length(x) > 0))
    if (has_any) {
        
        par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 3, 0))
        
        display_names <- list(
            outcomes  = c(masld_numeric = "MASLD", slf_binary = "SLF"),
            mediators = c(METS_IR = "METS-IR", HOMA_IR = "HOMA-IR")
        )
        
        color_palette <- list(
            predictor      = "#A9CCE3",
            outcome_masld  = "#F5B7B1",
            outcome_slf    = "#D7BDE2",
            mediator_mets  = "#FAD7A0",
            mediator_homa  = "#A9DFBF"
        )
        
        plot_counter <- 1
        plot_letters <- c("A", "B", "C", "D")
        
        for (outcome_var in names(mediation_results)) {
            for (mediator_var in names(mediation_results[[outcome_var]])) {
                
                mediation_bundle <- mediation_results[[outcome_var]][[mediator_var]]
                if (is.null(mediation_bundle)) next
                
                current_colors <- list(
                    predictor = color_palette$predictor,
                    mediator  = ifelse(mediator_var == "METS_IR", color_palette$mediator_mets, color_palette$mediator_homa),
                    outcome   = ifelse(outcome_var == "masld_numeric", color_palette$outcome_masld, color_palette$outcome_slf)
                )
                
                plot_publication_mediation(
                    mediation_bundle = mediation_bundle,
                    predictor_var    = "ALB_NLR",
                    predictor_label  = "Albumin/NLR",
                    mediator_var     = mediator_var,
                    mediator_label   = display_names$mediators[mediator_var],
                    outcome_label    = display_names$outcomes[outcome_var],
                    title_letter     = plot_letters[plot_counter],
                    color_palette    = current_colors
                )
                
                plot_counter <- plot_counter + 1
            }
        }
        
        mtext("Mediation Analyses of Albumin/NLR on Liver Outcomes", outer = TRUE, cex = 1.5, font = 2)
        
        par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0), family = "sans")
        cat("--- Plotting complete. Check the Plots pane for the diagrams. ---\n")
        
    } else {
        cat("--- No valid mediation results were available to plot. ---\n")
    }
    
    ```
    
- **Relationship CURI and Potential Mediators mediaton(bảng)**
    
    
    ```r
    # ==============================================================================
    # END-TO-END (NO SPLIT):
    # Outcomes: Alb/NLR (continuous), MASLD (binary), SLF (binary within MASLD)
    # Predictors: HOMA_IR + METS_IR
    # Models: M1 unadjusted; M2 + sociodemographics; M3 + comorbidities
    # ==============================================================================
    suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(stringr)
      library(gtsummary)
      library(gt)
    })
    
    if (!exists("MASLD")) stop("Object `MASLD` not found.")
    data <- as.data.frame(MASLD)
    
    ensure_cols <- function(df, cols){
      miss <- setdiff(cols, names(df))
      if (length(miss)) for (m in miss) df[[m]] <- NA
      df
    }
    
    calculate_alcohol_consumption <- function(data) {
      data <- ensure_cols(data, c("ALQ121","ALQ130","RIAGENDR"))
      data %>%
        mutate(
          days_per_week_alcohol = case_when(
            is.na(ALQ121) ~ NA_real_,
            ALQ121 == 0 ~ 0,
            ALQ121 == 1 ~ 7,
            ALQ121 == 2 ~ 6,
            ALQ121 == 3 ~ 3.5,
            ALQ121 == 4 ~ 2,
            ALQ121 == 5 ~ 1,
            ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
            ALQ121 == 7 ~ 1 / (30.4375 / 7),
            ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
            ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
            ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
            TRUE ~ NA_real_
          ),
          ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, as.numeric(ALQ130)),
          weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
        )
    }
    
    apply_exclusion_criteria <- function(data){
      cat("--- START: APPLYING EXCLUSION CRITERIA ---\n")
      n0 <- nrow(data)
    
      data <- ensure_cols(data, c(
        "RHD143","LBXHBC","LBDHBG","LBXHCR",
        "LUXCAPM","LUXSMED",
        "ALQ121","ALQ130","RIAGENDR","RIDAGEYR",
        "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"
      ))
    
      rows_before <- nrow(data)
      hep_vars <- c("RHD143","LBXHBC","LBDHBG","LBXHCR")
      data <- data %>% filter(rowSums(across(all_of(hep_vars), ~ .x == 1), na.rm = TRUE) == 0)
      cat("1) Excluded viral hepatitis:", rows_before - nrow(data), "\n")
    
      rows_before <- nrow(data)
      data <- data %>%
        filter(!is.na(LUXCAPM), !is.na(LUXSMED)) %>%
        filter(!is.na(LBDNENO), !is.na(LBDLYMNO)) %>%
        filter(!is.na(LBXNEPCT), !is.na(LBXLYPCT)) %>%
        filter(!(is.na(LBXSAL) & is.na(LBDSALSI)))
      cat("2) Excluded missing CAP/LSM/blood:", rows_before - nrow(data), "\n")
    
      data <- calculate_alcohol_consumption(data)
      rows_before <- nrow(data)
      data <- data %>%
        filter(
          is.na(weekly_alcohol_grams) |
            (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
            (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
        )
      cat("3) Excluded heavy alcohol:", rows_before - nrow(data), "\n")
    
      rows_before <- nrow(data)
      data <- data %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
      cat("4) Excluded age <18:", rows_before - nrow(data), "\n")
    
      cat("--- EXCLUSION COMPLETE. Removed:", n0 - nrow(data), "Remaining:", nrow(data), "---\n\n")
      data
    }
    
    create_all_features <- function(data){
      cat("--- START: FEATURE ENGINEERING ---\n")
    
      data <- ensure_cols(data, c(
        "LUXCAPM","LUXSMED",
        "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI",
        "LBDGLUSI","LBDINSI","LBXGLU","LBXTR","BMXBMI","LBDHDD",
        "BMXWAIST","DIQ010","DIQ050","DIQ070","LBXGH","LBXSGL",
        "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3","BPQ040A","BPQ020",
        "LBDTRSI","LBDHDDSI","BPQ090D",
        "RIDAGEYR","RIAGENDR","RIDRETH1","RIDRETH3","DMDEDUC2","SMQ020","SMQ040",
        "MCQ220","KIQ022","KIQ025","LBDSCRSI","LBXSCR",
        "MCQ160B","MCQ160C","MCQ160E","MCQ160F"
      ))
    
      out <- data %>%
        mutate(
          albumin_g_dl = case_when(
            !is.na(LBXSAL) ~ as.numeric(LBXSAL),
            is.na(LBXSAL) & !is.na(LBDSALSI) ~ as.numeric(LBDSALSI) / 10,
            TRUE ~ NA_real_
          ),
          neutrophil_count = as.numeric(LBDNENO),
          lymphocyte_count = as.numeric(LBDLYMNO),
          NLR = if_else(!is.na(neutrophil_count) & !is.na(lymphocyte_count) & lymphocyte_count > 0,
                       neutrophil_count / lymphocyte_count, NA_real_),
          Alb_by_NLR = if_else(!is.na(albumin_g_dl) & !is.na(NLR) & NLR > 0, albumin_g_dl / NLR, NA_real_),
    
          HOMA_IR = if_else(!is.na(LBDGLUSI) & !is.na(LBDINSI),
                           (as.numeric(LBDGLUSI) * as.numeric(LBDINSI)) / 22.5, NA_real_),
    
          METS_IR = if_else(!is.na(LBXGLU) & !is.na(LBXTR) & !is.na(BMXBMI) & !is.na(LBDHDD) &
                              as.numeric(LBDHDD) > 0 & log(as.numeric(LBDHDD)) != 0,
                            (log(2 * as.numeric(LBXGLU) + as.numeric(LBXTR)) * as.numeric(BMXBMI)) / log(as.numeric(LBDHDD)),
                            NA_real_)
        ) %>%
        mutate(
          has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
    
          risk_bmi_waist = if_else(
            is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
            NA,
            (BMXBMI >= 25) |
              (RIAGENDR == 1 & BMXWAIST >= 94) |
              (RIAGENDR == 2 & BMXWAIST >= 80)
          ),
    
          risk_glucose_diabetes = if_else(
            is.na(LBXGLU) & is.na(LBXGH) &
              (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
              (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
              (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
            NA,
            (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
              (if_else(is.na(LBXGH),  FALSE, LBXGH >= 5.7)) |
              (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
              (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
              (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
          ),
    
          any_bp_measurement_high =
            (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
            (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
            (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
    
          risk_blood_pressure = if_else(
            is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
            NA,
            any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
          ),
    
          risk_triglycerides = if_else(
            is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
            NA,
            (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          ),
    
          risk_low_hdl = if_else(
            is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))),
            NA,
            ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
               (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          )
        ) %>%
        rowwise() %>%
        mutate(num_cardiometabolic_risks = sum(c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure,
                                                risk_triglycerides, risk_low_hdl), na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(
          has_one_plus_cardiometabolic_risk = if_else(
            is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
              is.na(risk_triglycerides) & is.na(risk_low_hdl),
            NA,
            num_cardiometabolic_risks > 0
          ),
          masld_group = case_when(
            has_hepatic_steatosis == FALSE ~ "non-MASLD",
            has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
            has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
            TRUE ~ NA_character_
          ),
          masld_binary = if_else(masld_group == "MASLD", 1, 0),
          slf_group = case_when(
            masld_group != "MASLD" ~ NA_character_,
            is.na(LUXSMED) ~ NA_character_,
            LUXSMED >= 8.0 ~ "SLF",
            LUXSMED < 8.0  ~ "non-SLF",
            TRUE ~ NA_character_
          ),
          slf_binary = if_else(slf_group == "SLF", 1, 0)
        ) %>%
        mutate(
          age = RIDAGEYR,
          gender = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
          race = factor(
            RIDRETH3,
            levels = c(1, 2, 3, 4, 6),
            labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian")
          ),
          education = factor(case_when(
            DMDEDUC2 == 1 ~ "Less than 9th grade",
            DMDEDUC2 == 2 ~ "9-11th grade",
            DMDEDUC2 == 3 ~ "High school graduate/GED",
            DMDEDUC2 == 4 ~ "Some college or AA degree",
            DMDEDUC2 == 5 ~ "College graduate or above",
            TRUE ~ NA_character_
          ), levels = c("Less than 9th grade", "9-11th grade", "High school graduate/GED",
                       "Some college or AA degree", "College graduate or above")),
          smoking_status = factor(case_when(
            SMQ020 == 2 ~ "Non-smoker",
            SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
            SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
            TRUE ~ NA_character_
          ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
          obesity = factor(if_else(!is.na(BMXBMI) & BMXBMI >= 30, "Yes", "No", missing = "No"),
                           levels = c("No", "Yes")),
          cancer = factor(if_else(MCQ220 == 1, "Yes", "No", missing = "No"), levels = c("No", "Yes")),
          diabetes_mellitus = factor(
            if_else(
              DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                (!is.na(LBXGH) & LBXGH >= 6.5) |
                (!is.na(LBXSGL) & LBXSGL >= 126),
              "Yes", "No", missing = "No"
            ),
            levels = c("No", "Yes")
          ),
          hypertension = factor(
            if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"),
            levels = c("No", "Yes")
          ),
          eGFR = {
            creatinine_mg_dl <- case_when(
              !is.na(LBXSCR) ~ as.numeric(LBXSCR),
              is.na(LBXSCR) & !is.na(LBDSCRSI) ~ as.numeric(LBDSCRSI) / 88.4,
              TRUE ~ NA_real_
            )
            if_else(!is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                    175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                      if_else(RIAGENDR == 2, 0.742, 1) * if_else(RIDRETH1 == 4, 1.212, 1),
                    NA_real_)
          },
          chronic_kidney_disease = factor(
            if_else(KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60), "Yes", "No", missing = "No"),
            levels = c("No", "Yes")
          ),
          Cardiovascular_Disease = factor(
            if_else(MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160F == 1, "Yes", "No", missing = "No"),
            levels = c("No", "Yes")
          )
        )
    
      cat("--- FEATURE ENGINEERING COMPLETE. ---\n\n")
      out
    }
    
    # RUN PREP
    data_cleaned <- apply_exclusion_criteria(data)
    data_processed <- create_all_features(data_cleaned)
    
    # -------------------------
    # 4) MODELING + TABLES
    # -------------------------
    preds <- c("HOMA_IR","METS_IR")
    cov_soc <- c("age","gender","race","education","smoking_status")
    cov_comorb <- c(cov_soc, "obesity","diabetes_mellitus","hypertension","cancer","chronic_kidney_disease","Cardiovascular_Disease")
    
    fit_cc <- function(df, formula, outcome_var, family = NULL){
      vars <- all.vars(formula)
      dcc <- df[, vars, drop = FALSE]
      dcc <- stats::na.omit(dcc)
      if (nrow(dcc) < (length(vars) + 10)) return(NULL)
      if (!is.null(family) && dplyr::n_distinct(dcc[[outcome_var]]) < 2) return(NULL)
      fcols <- names(dcc)[vapply(dcc, is.factor, logical(1))]
      for (cc in fcols) if (dplyr::n_distinct(dcc[[cc]]) < 2) return(NULL)
      if (is.null(family)) stats::lm(formula, data = dcc) else stats::glm(formula, data = dcc, family = family)
    }
    
    make_block <- function(df, outcome_var, type = c("linear","logistic"), header_label){
      type <- match.arg(type)
    
      if (outcome_var == "slf_binary") df <- df %>% filter(masld_group == "MASLD")
    
      p_ok <- intersect(preds, names(df))
      c2 <- intersect(cov_soc, names(df))
      c3 <- intersect(cov_comorb, names(df))
    
      f1 <- as.formula(paste0(outcome_var, " ~ ", paste(p_ok, collapse=" + ")))
      f2 <- as.formula(paste0(outcome_var, " ~ ", paste(c(p_ok, c2), collapse=" + ")))
      f3 <- as.formula(paste0(outcome_var, " ~ ", paste(c(p_ok, c3), collapse=" + ")))
    
      if (type == "linear") {
        m1 <- fit_cc(df, f1, outcome_var, family = NULL)
        m2 <- fit_cc(df, f2, outcome_var, family = NULL)
        m3 <- fit_cc(df, f3, outcome_var, family = NULL)
    
        t1 <- if (!is.null(m1)) tbl_regression(m1, include = tidyselect::all_of(p_ok)) %>% bold_p(t=0.05) else NULL
        t2 <- if (!is.null(m2)) tbl_regression(m2, include = tidyselect::all_of(p_ok)) %>% bold_p(t=0.05) else NULL
        t3 <- if (!is.null(m3)) tbl_regression(m3, include = tidyselect::all_of(p_ok)) %>% bold_p(t=0.05) else NULL
      } else {
        m1 <- fit_cc(df, f1, outcome_var, family = binomial(link="logit"))
        m2 <- fit_cc(df, f2, outcome_var, family = binomial(link="logit"))
        m3 <- fit_cc(df, f3, outcome_var, family = binomial(link="logit"))
    
        t1 <- if (!is.null(m1)) tbl_regression(m1, exponentiate = TRUE, include = tidyselect::all_of(p_ok)) %>% bold_p(t=0.05) else NULL
        t2 <- if (!is.null(m2)) tbl_regression(m2, exponentiate = TRUE, include = tidyselect::all_of(p_ok)) %>% bold_p(t=0.05) else NULL
        t3 <- if (!is.null(m3)) tbl_regression(m3, exponentiate = TRUE, include = tidyselect::all_of(p_ok)) %>% bold_p(t=0.05) else NULL
      }
    
      tbls <- list(t1,t2,t3)
      spanners <- c("**Model 1**","**Model 2**","**Model 3**")
      keep <- !sapply(tbls, is.null)
      if (!any(keep)) return(NULL)
    
      merged <- tbl_merge(tbls = tbls[keep], tab_spanner = spanners[keep])
      list(tbl = merged, header = header_label)
    }
    
    cat("--- START: MODELS (NO SPLIT) ---\n")
    
    blk_albnlr <- make_block(data_processed, "Alb_by_NLR", "linear",  "**Outcome: Alb/NLR**")
    blk_masld  <- make_block(data_processed, "masld_binary", "logistic","**Outcome: MASLD**")
    blk_slf    <- make_block(data_processed, "slf_binary", "logistic","**Outcome: SLF**")
    
    blocks <- list(blk_albnlr, blk_masld, blk_slf)
    blocks <- blocks[!sapply(blocks, is.null)]
    if (length(blocks) == 0) stop("CRITICAL: No tables produced (insufficient complete cases or no outcome variation).")
    
    final_tbl <- tbl_stack(
      tbls = lapply(blocks, `[[`, "tbl"),
      group_header = sapply(blocks, `[[`, "header"),
      quiet = TRUE
    ) %>%
      modify_caption("**Table. Association of HOMA-IR and METS-IR with Alb/NLR, MASLD, and SLF (no split)**")
    
    final_gt <- final_tbl %>%
      as_gt() %>%
      tab_source_note(source_note = md("Abbreviations: CI = Confidence interval; OR = Odds ratio; Alb/NLR = albumin_g_dl / NLR.")) %>%
      tab_source_note(source_note = md("**Model 1**: Unadjusted. **Model 2**: + sociodemographics (Age, Gender, Race, Education, Smoking status). **Model 3**: Model 2 + comorbidities (Obesity, Diabetes mellitus, Hypertension, Cancer, Chronic kidney disease, Cardiovascular disease)."))
    
    print(final_gt)
    
    cat("\n--- COMPLETE ---\n")
    
    ```
    
- **AUC**
    
    ```r
    # ==============================================================================
    # END-TO-END: ROC + AUC/CI + Youden cutoff + LR+/LR- + DeLong test (NO split)
    # + ADD: Combined-score logistic models (e.g., Alb/NLR + HSI) with:
    #   - predicted probability columns
    #   - coefficient table (Beta/SE/p/OR/95%CI)
    #   - equation strings (logit(p)=...)
    # Predictors ONLY:
    #   Alb/NLR + HSI, HSI, FIB-4, Alb/NLR, SII, SIRI, AISI, NAR, PNI
    # Outcomes:
    #   MASLD (whole eligible), SLF (within MASLD only)
    # Remove: CURI, hs-CRP, uric acid
    # NHANES blood vars (primary):
    #   LBXSAL (Alb g/dL), LBXLCC (Lymph 10^3/uL), LBXNEU (Neut 10^3/uL),
    #   LBXPLT or LBXPLTSI (Platelet 10^3/uL), LBXMOC (Mono 10^3/uL),
    #   LBXNEPCT/LBXLYPCT (% fallback), LBDNENO/LBDLYMNO (count fallback), LBDSALSI (Alb SI)
    # FIX:
    #   SII = (Neut_count_k * Plt_count_k / Lymph_count_k) * 1000
    #   Plt_count_k = coalesce(LBXPLT, LBXPLTSI)
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(pROC)
        library(ggplot2)
        library(patchwork)
        library(RColorBrewer)
        library(stringr)
        library(grid)
        library(broom)
    })
    
    # ----------------------------- 0) HELPERS -----------------------------
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    ensure_cols <- function(df, cols, default = NA) {
        miss <- setdiff(cols, names(df))
        if (length(miss)) for (m in miss) df[[m]] <- default
        df
    }
    
    safe_div <- function(a, b) {
        ifelse(is.na(a) | is.na(b) | !is.finite(a) | !is.finite(b) | b == 0, NA_real_, a / b)
    }
    
    get_predictor_name <- function(p_var, outcome) {
        ifelse(grepl("\\+", p_var),
               gsub("\\+", "_", paste0(outcome, "_", p_var)),
               p_var)
    }
    
    # ----------------------------- 1) ALCOHOL -----------------------------
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_cols(data, c("ALQ121","ALQ130","RIAGENDR"), default = NA)
        data %>%
            mutate(
                ALQ121 = safe_num(ALQ121),
                ALQ130 = safe_num(ALQ130),
                RIAGENDR = safe_num(RIAGENDR),
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2,3)) / (30.4375/7),
                    ALQ121 == 7 ~ 1 / (30.4375/7),
                    ALQ121 == 8 ~ mean(c(7,11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3,6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1,2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777,999) | is.na(ALQ130), NA_real_, ALQ130),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # ----------------------------- 2) EXCLUSION / CLEAN -----------------------------
    clean_initial_data <- function(data) {
        cat("--- START: EXCLUSION / CLEANING ---\n")
        initial_rows_total <- nrow(data)
        
        data <- data[, colSums(!is.na(data)) > 0, drop = FALSE]
        
        data <- ensure_cols(
            data,
            c(
                "RIDAGEYR","RIAGENDR",
                "LUXSMED","LUXCAPM",
                "ALQ121","ALQ130",
                "RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV",
                # Albumin
                "LBXSAL","LBDSALSI",
                # Primary counts/% for indices
                "LBXLCC","LBXNEU","LBXNEPCT",
                # Platelets
                "LBXPLT","LBXPLTSI",
                "LBXMOC",
                # Fallbacks
                "LBXLYPCT","LBDLYMNO","LBDNENO"
            ),
            default = NA
        )
        
        cat("1) Viral hepatitis exclusion...\n")
        hep_vars <- c("RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV")
        hep_vars_exist <- intersect(hep_vars, names(data))
        before_n <- nrow(data)
        if (length(hep_vars_exist) > 0) {
            data <- data %>%
                mutate(across(all_of(hep_vars_exist), safe_num)) %>%
                filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
        }
        cat("-> Removed:", before_n - nrow(data), "\n\n")
        
        cat("2) Require FibroScan + Albumin + Neut/Lymph (counts OR %)...\n")
        before_n <- nrow(data)
        
        data <- data %>% drop_na(LUXSMED, LUXCAPM)
        
        data <- data %>%
            mutate(
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                Albumin_g_dL_tmp = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                
                LBXNEU   = safe_num(LBXNEU),
                LBXLCC   = safe_num(LBXLCC),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXLYPCT = safe_num(LBXLYPCT),
                
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                
                has_counts_primary  = !is.na(LBXNEU) & !is.na(LBXLCC),
                has_counts_fallback = !is.na(LBDNENO) & !is.na(LBDLYMNO),
                has_pct             = !is.na(LBXNEPCT) & !is.na(LBXLYPCT)
            ) %>%
            filter(!is.na(Albumin_g_dL_tmp)) %>%
            filter(has_counts_primary | has_counts_fallback | has_pct)
        
        cat("-> Removed:", before_n - nrow(data), "\n\n")
        
        cat("3) Alcohol threshold exclusion...\n")
        data <- calculate_alcohol_consumption(data)
        before_n <- nrow(data)
        data <- data %>%
            mutate(RIAGENDR = safe_num(RIAGENDR)) %>%
            filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        cat("-> Removed:", before_n - nrow(data), "\n\n")
        
        cat("4) Age >= 18...\n")
        before_n <- nrow(data)
        data <- data %>%
            mutate(RIDAGEYR = safe_num(RIDAGEYR)) %>%
            filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        cat("-> Removed:", before_n - nrow(data), "\n\n")
        
        cat("--- DONE. Total removed:", initial_rows_total - nrow(data),
            "| Remaining:", nrow(data), "---\n\n")
        data
    }
    
    # ----------------------------- 3) FEATURE ENGINEERING -----------------------------
    create_features <- function(data) {
        cat("--- START: FEATURE ENGINEERING ---\n")
        
        cols_need <- c(
            "BMXBMI","BMXWAIST","BMXHIP",
            "LBXGLU","LBXGH","LBXSGL",
            "DIQ010","DIQ050","DIQ070",
            "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
            "BPQ040A","BPQ090D","BPQ020",
            "LBDTRSI","LBDHDDSI","LBDTCSI",
            "LBXSATSI","LBXSASSI","LBXPLTSI","LBXPLT",
            "LBXALT","LBXAST"
        )
        data <- ensure_cols(data, cols_need, default = NA)
        
        data2 <- data %>%
            mutate(
                RIAGENDR = safe_num(RIAGENDR),
                RIDAGEYR = safe_num(RIDAGEYR),
                LUXCAPM  = safe_num(LUXCAPM),
                LUXSMED  = safe_num(LUXSMED),
                
                BMXBMI   = safe_num(BMXBMI),
                BMXWAIST = safe_num(BMXWAIST),
                BMXHIP   = safe_num(BMXHIP),
                LBXGLU   = safe_num(LBXGLU),
                LBXGH    = safe_num(LBXGH),
                LBXSGL   = safe_num(LBXSGL),
                
                BPXOSY1  = safe_num(BPXOSY1),
                BPXODI1  = safe_num(BPXODI1),
                BPXOSY2  = safe_num(BPXOSY2),
                BPXODI2  = safe_num(BPXODI2),
                BPXOSY3  = safe_num(BPXOSY3),
                BPXODI3  = safe_num(BPXODI3),
                
                BPQ040A  = safe_num(BPQ040A),
                BPQ090D  = safe_num(BPQ090D),
                BPQ020   = safe_num(BPQ020),
                
                LBDTRSI  = safe_num(LBDTRSI),
                LBDHDDSI = safe_num(LBDHDDSI),
                LBDTCSI  = safe_num(LBDTCSI),
                
                LBXSATSI = safe_num(LBXSATSI),
                LBXSASSI = safe_num(LBXSASSI),
                LBXALT   = safe_num(LBXALT),
                LBXAST   = safe_num(LBXAST),
                
                # platelets: both
                LBXPLTSI = safe_num(LBXPLTSI),
                LBXPLT   = safe_num(LBXPLT),
                
                # albumin + CBC
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                LBXNEU   = safe_num(LBXNEU),
                LBXLCC   = safe_num(LBXLCC),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXMOC   = safe_num(LBXMOC),
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBXLYPCT = safe_num(LBXLYPCT)
            ) %>%
            mutate(
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                
                risk_bmi_waist = if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                    NA,
                    (BMXBMI >= 25) |
                        (RIAGENDR == 1 & BMXWAIST >= 94) |
                        (RIAGENDR == 2 & BMXWAIST >= 80)
                ),
                
                risk_glucose_diabetes = if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7,9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7,9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7,9)),
                    NA,
                    (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
                        (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                
                any_bp_high =
                    (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                
                risk_blood_pressure = if_else(
                    is.na(any_bp_high) & (is.na(BPQ040A) | BPQ040A %in% c(7,9)),
                    NA,
                    any_bp_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                
                risk_triglycerides = if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9)),
                    NA,
                    (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                
                risk_low_hdl = if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9))),
                    NA,
                    ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                )
            ) %>%
            rowwise() %>%
            mutate(
                num_cardiometabolic_risks = sum(
                    c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
                    na.rm = TRUE
                )
            ) %>%
            ungroup() %>%
            mutate(
                has_one_plus_cardiometabolic_risk = if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA,
                    num_cardiometabolic_risks > 0
                ),
                
                is_light_drinker_final = case_when(
                    is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
                    (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                
                masld_group = case_when(
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
                    TRUE ~ NA_character_
                ),
                masld_binary = if_else(masld_group == "MASLD", 1L, 0L, missing = NA_integer_),
                
                slf_group = case_when(
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ "SLF",
                    masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED < 8.0  ~ "non-SLF",
                    TRUE ~ NA_character_
                ),
                slf_binary = if_else(slf_group == "SLF", 1L, 0L, missing = NA_integer_)
            )
        
        # ===================== BLOOD DERIVED INDICES =====================
        data2 <- data2 %>%
            mutate(
                Albumin_g_dL = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                
                # counts (10^3/uL)
                Neut_count_k  = dplyr::coalesce(LBXNEU,  LBDNENO),
                Lymph_count_k = dplyr::coalesce(LBXLCC,  LBDLYMNO),
                Mono_count_k  = LBXMOC,
                
                # platelet count (10^3/uL): FIX fallback
                Plt_count_k   = dplyr::coalesce(LBXPLT, LBXPLTSI),
                
                Neut_percent  = LBXNEPCT,
                Lymph_percent = LBXLYPCT,
                
                # NLR: counts preferred; else % ratio
                NLR = case_when(
                    !is.na(Neut_count_k) & !is.na(Lymph_count_k) & Lymph_count_k > 0 ~ Neut_count_k / Lymph_count_k,
                    (is.na(Neut_count_k) | is.na(Lymph_count_k)) &
                        !is.na(Neut_percent) & !is.na(Lymph_percent) & Lymph_percent > 0 ~ Neut_percent / Lymph_percent,
                    TRUE ~ NA_real_
                ),
                
                Alb_NLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                                  Albumin_g_dL / NLR, NA_real_),
                
                # PNI = 10*Alb + 0.005*(TLC*1000)
                PNI = if_else(!is.na(Albumin_g_dL) & !is.na(Lymph_count_k),
                              10 * Albumin_g_dL + 0.005 * (Lymph_count_k * 1000),
                              NA_real_),
                
                # NAR = (ANC/Alb)*1000
                NAR = if_else(!is.na(Neut_count_k) & !is.na(Albumin_g_dL) & Albumin_g_dL > 0,
                              (Neut_count_k / Albumin_g_dL) * 1000,
                              NA_real_),
                
                # SII = (N*P/L)*1000
                SII = if_else(!is.na(Neut_count_k) & !is.na(Plt_count_k) & !is.na(Lymph_count_k) & Lymph_count_k > 0,
                              (Neut_count_k * Plt_count_k / Lymph_count_k) * 1000,
                              NA_real_),
                
                # SIRI = (N*M/L)*1000
                SIRI = if_else(!is.na(Neut_count_k) & !is.na(Mono_count_k) & !is.na(Lymph_count_k) & Lymph_count_k > 0,
                               (Neut_count_k * Mono_count_k / Lymph_count_k) * 1000,
                               NA_real_),
                
                # AISI = (N*P*M/L)*1e6
                AISI = if_else(!is.na(Neut_count_k) & !is.na(Plt_count_k) & !is.na(Mono_count_k) &
                                   !is.na(Lymph_count_k) & Lymph_count_k > 0,
                               (Neut_count_k * Plt_count_k * Mono_count_k / Lymph_count_k) * 1e6,
                               NA_real_)
            )
        
        # ===================== HSI + FIB-4 =====================
        data2 <- data2 %>%
            mutate(
                gender = factor(RIAGENDR, labels = c("Male","Female")),
                
                diabetes_mellitus = factor(
                    if_else(
                        DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                            (!is.na(LBXGH)  & LBXGH  >= 6.5) |
                            (!is.na(LBXSGL) & LBXSGL >= 126),
                        "Yes","No", missing = "No"
                    ),
                    levels = c("No","Yes")
                ),
                
                ALT_u = dplyr::coalesce(LBXSATSI, LBXALT),
                AST_u = dplyr::coalesce(LBXSASSI, LBXAST),
                
                PLT_k_for_fib4 = dplyr::coalesce(LBXPLT, LBXPLTSI),
                
                FIB4 = if_else(
                    !is.na(RIDAGEYR) & !is.na(AST_u) & !is.na(ALT_u) & !is.na(PLT_k_for_fib4) &
                        PLT_k_for_fib4 > 0 & ALT_u > 0,
                    (RIDAGEYR * AST_u) / (PLT_k_for_fib4 * sqrt(ALT_u)),
                    NA_real_
                ),
                
                HSI = if_else(
                    !is.na(ALT_u) & !is.na(AST_u) & AST_u != 0 & !is.na(BMXBMI),
                    8 * (ALT_u / AST_u) + BMXBMI +
                        (2 * if_else(diabetes_mellitus == "Yes", 1, 0, missing = 0)) +
                        (2 * if_else(gender == "Female", 1, 0, missing = 0)),
                    NA_real_
                )
            )
        
        cat("--- DONE: FEATURE ENGINEERING ---\n")
        cat("Non-missing counts: SII =", sum(!is.na(data2$SII) & is.finite(data2$SII)),
            "| Neut_count_k =", sum(!is.na(data2$Neut_count_k)),
            "| Plt_count_k =", sum(!is.na(data2$Plt_count_k)),
            "| Lymph_count_k =", sum(!is.na(data2$Lymph_count_k)), "\n\n")
        data2
    }
    
    # ----------------------------- 4) COMBINED SCORE MODELS + EQUATIONS + COEFS -----------------------------
    extract_model_stats <- function(model, outcome_name, model_name) {
        broom::tidy(model) %>%
            mutate(
                Outcome = outcome_name,
                Model   = model_name,
                OR      = exp(estimate),
                CI_Lower = exp(estimate - 1.96 * std.error),
                CI_Upper = exp(estimate + 1.96 * std.error)
            ) %>%
            select(Outcome, Model, Term = term, Beta = estimate, SE = std.error, P_value = p.value,
                   OR, CI_Lower, CI_Upper)
    }
    
    get_equation_string <- function(model, digits = 6) {
        cf <- coef(model)
        if (length(cf) == 0) return(NA_character_)
        intercept <- unname(cf[1])
        preds <- cf[-1]
        eq <- sprintf(paste0("logit(p) = %.", digits, "f"), intercept)
        if (length(preds)) {
            for (i in seq_along(preds)) {
                b <- unname(preds[i])
                nm <- names(preds)[i]
                eq <- paste0(eq, sprintf(paste0(" + (%.", digits, "f * %s)"), b, nm))
            }
        }
        eq
    }
    
    fit_combined_models_and_write <- function(df, outcomes, combos,
                                              min_n = 40,
                                              out_eq_csv = "CombinedModels_Equations.csv",
                                              out_coef_csv = "CombinedModels_Coefficients.csv") {
        model_stats_list <- list()
        model_equations_list <- list()
        
        for (outcome in outcomes) {
            for (cmb in combos) {
                parts <- strsplit(cmb, "\\+")[[1]] %>% trimws()
                if (length(parts) != 2) next
                x1 <- parts[1]; x2 <- parts[2]
                if (!all(c(outcome, x1, x2) %in% names(df))) next
                
                # SLF models: fit/predict within MASLD only, keep others NA
                if (outcome == "slf_binary") {
                    fit_base <- df %>% filter(masld_group == "MASLD")
                    pred_rows <- which(df$masld_group == "MASLD")
                } else {
                    fit_base <- df
                    pred_rows <- seq_len(nrow(df))
                }
                
                d_fit <- fit_base %>%
                    filter(!is.na(.data[[outcome]]), .data[[outcome]] %in% c(0,1)) %>%
                    filter(!is.na(.data[[x1]]), is.finite(.data[[x1]])) %>%
                    filter(!is.na(.data[[x2]]), is.finite(.data[[x2]]))
                
                if (nrow(d_fit) < min_n || dplyr::n_distinct(d_fit[[outcome]]) < 2) next
                
                m <- tryCatch(
                    glm(as.formula(paste(outcome, "~", x1, "+", x2)),
                        data = d_fit, family = binomial()),
                    error = function(e) NULL
                )
                if (is.null(m)) next
                
                # 1) predicted prob column name: outcome_Alb_NLR+HSI -> outcome_Alb_NLR_HSI
                pred_col <- get_predictor_name(cmb, outcome)
                if (!pred_col %in% names(df)) df[[pred_col]] <- NA_real_
                
                pr <- tryCatch(predict(m, newdata = df[pred_rows, , drop = FALSE], type = "response"),
                               error = function(e) rep(NA_real_, length(pred_rows)))
                df[[pred_col]][pred_rows] <- pr
                
                # 2) store coefficient table
                model_name_full <- paste0(str_replace_all(cmb, "\\+", " + "))
                model_stats_list[[paste(outcome, cmb, sep = "__")]] <- extract_model_stats(m, outcome, model_name_full)
                
                # 3) store equation text
                eq <- get_equation_string(m, digits = 6)
                full_sentence <- sprintf("For the prediction of %s, the model was %s for %s.",
                                         outcome, eq, model_name_full)
                
                model_equations_list[[paste(outcome, cmb, sep = "__")]] <- data.frame(
                    Outcome = outcome,
                    Model   = model_name_full,
                    Equation = eq,
                    Equation_Text = full_sentence,
                    stringsAsFactors = FALSE
                )
            }
        }
        
        equations_df <- dplyr::bind_rows(model_equations_list)
        coeffs_df <- dplyr::bind_rows(model_stats_list)
        
        # write outputs (safe even if empty)
        write.csv(equations_df, out_eq_csv, row.names = FALSE)
        write.csv(coeffs_df, out_coef_csv, row.names = FALSE)
        
        cat("--- Combined models outputs ---\n")
        cat("Saved:", out_eq_csv, "\n")
        cat("Saved:", out_coef_csv, "\n\n")
        
        if (nrow(equations_df) > 0) {
            cat(">>> EQUATIONS (Copy paste):\n")
            for (i in seq_len(nrow(equations_df))) {
                cat(paste0("- ", equations_df$Equation_Text[i], "\n\n"))
            }
        } else {
            cat(">>> EQUATIONS: none generated (insufficient data after filters)\n\n")
        }
        
        list(data = df, equations = equations_df, coefficients = coeffs_df)
    }
    
    # ----------------------------- 5) ROC ANALYSIS + DeLong -----------------------------
    analyze_roc_with_delong <- function(data, outcome, predictors, ref_predictor) {
        results <- list()
        
        for (p_var in predictors) {
            current_p <- get_predictor_name(p_var, outcome)
            if (!current_p %in% names(data) || !outcome %in% names(data)) next
            
            valid_idx <- !is.na(data[[outcome]]) & data[[outcome]] %in% c(0,1) &
                !is.na(data[[current_p]]) & is.finite(data[[current_p]])
            if (sum(valid_idx) < 20) next
            
            d <- data[valid_idx, , drop = FALSE]
            if (length(unique(d[[outcome]])) < 2 || length(unique(d[[current_p]])) < 2) next
            
            roc_obj <- tryCatch(
                roc(d[[outcome]], d[[current_p]], quiet = TRUE, levels = c(0,1), direction = "<"),
                error = function(e) NULL
            )
            if (is.null(roc_obj)) next
            
            auc_val <- as.numeric(auc(roc_obj))
            ci_val  <- tryCatch(ci.auc(roc_obj), error = function(e) c(NA, NA, NA))
            
            opt <- tryCatch(
                coords(roc_obj, "best",
                       ret = c("threshold","sensitivity","specificity"),
                       best.method = "youden"),
                error = function(e) NULL
            )
            if (is.null(opt) || nrow(opt) == 0) next
            
            sens <- as.numeric(opt[1,"sensitivity"])
            spec <- as.numeric(opt[1,"specificity"])
            lr_plus  <- sens / (1 - spec + .Machine$double.eps)
            lr_minus <- (1 - sens) / (spec + .Machine$double.eps)
            
            if (!is.na(auc_val) && auc_val < 0.5) {
                auc_val <- 1 - auc_val
                if (all(!is.na(ci_val))) ci_val <- 1 - rev(ci_val)
            }
            
            results[[p_var]] <- data.frame(
                Outcome = outcome,
                Predictor = p_var,
                AUC = round(auc_val, 3),
                CI  = if (any(is.na(ci_val))) NA_character_ else sprintf("%.3f–%.3f", ci_val[1], ci_val[3]),
                Cutoff = round(as.numeric(opt[1,"threshold"]), 3),
                Sensitivity = round(sens, 3),
                Specificity = round(spec, 3),
                LR_Plus  = round(lr_plus, 3),
                LR_Minus = round(lr_minus, 3),
                stringsAsFactors = FALSE
            )
        }
        
        if (length(results) == 0) return(data.frame())
        results_df <- dplyr::bind_rows(results)
        
        ref_name <- get_predictor_name(ref_predictor, outcome)
        if (!ref_name %in% names(data)) {
            results_df$p_value <- NA_real_
            return(results_df %>% arrange(desc(AUC)))
        }
        
        p_values <- sapply(predictors, function(p_var) {
            if (p_var == ref_predictor) return(NA_real_)
            cur_name <- get_predictor_name(p_var, outcome)
            if (!cur_name %in% names(data)) return(NA_real_)
            
            common_idx <- !is.na(data[[outcome]]) & data[[outcome]] %in% c(0,1) &
                !is.na(data[[ref_name]]) & is.finite(data[[ref_name]]) &
                !is.na(data[[cur_name]]) & is.finite(data[[cur_name]])
            
            if (sum(common_idx) < 20) return(NA_real_)
            d <- data[common_idx, , drop = FALSE]
            if (dplyr::n_distinct(d[[outcome]]) < 2) return(NA_real_)
            
            ref_roc <- tryCatch(roc(d[[outcome]], d[[ref_name]], quiet = TRUE), error = function(e) NULL)
            cur_roc <- tryCatch(roc(d[[outcome]], d[[cur_name]], quiet = TRUE), error = function(e) NULL)
            if (is.null(ref_roc) || is.null(cur_roc)) return(NA_real_)
            
            test_res <- tryCatch(roc.test(ref_roc, cur_roc, method = "delong"), error = function(e) NULL)
            if (is.null(test_res)) NA_real_ else test_res$p.value
        })
        
        results_df$p_value <- p_values[results_df$Predictor]
        results_df %>% arrange(desc(AUC))
    }
    
    plot_roc_curves_delong <- function(data, result_table, outcome_name, parameters, cohort_label, ref_label) {
        
        display_names <- c(
            "Alb_NLR+HSI" = "Alb/NLR + HSI",
            "HSI" = "HSI",
            "FIB4" = "FIB-4",
            "Alb_NLR" = "Alb/NLR",
            "SII" = "SII",
            "SIRI" = "SIRI",
            "AISI" = "AISI",
            "NAR" = "NAR",
            "PNI" = "PNI"
        )
        
        outcome_display_names <- c(
            "masld_binary" = "MASLD",
            "slf_binary"   = "Significant Liver Fibrosis (SLF)"
        )
        
        colour_palette <- brewer.pal(n = max(3, length(parameters)), "Dark2")[1:length(parameters)]
        names(colour_palette) <- parameters
        
        preds_for_outcome <- result_table$Predictor[result_table$Outcome == outcome_name]
        if (length(preds_for_outcome) == 0) return(NULL)
        
        plot_df_list <- lapply(preds_for_outcome, function(pred) {
            cur <- get_predictor_name(pred, outcome_name)
            if (!cur %in% names(data)) return(NULL)
            
            idx <- !is.na(data[[outcome_name]]) & data[[outcome_name]] %in% c(0,1) &
                !is.na(data[[cur]]) & is.finite(data[[cur]])
            if (sum(idx) < 20) return(NULL)
            
            roc_obj <- roc(data[[outcome_name]][idx], data[[cur]][idx], quiet = TRUE, levels = c(0,1), direction = "<")
            if (as.numeric(auc(roc_obj)) < 0.5) {
                roc_obj <- roc(roc_obj$response, -roc_obj$predictor, quiet = TRUE, levels = c(0,1), direction = "<")
            }
            
            data.frame(
                fpr = 1 - roc_obj$specificities,
                tpr = roc_obj$sensitivities,
                predictor = factor(pred, levels = parameters),
                stringsAsFactors = FALSE
            )
        })
        
        plot_df <- do.call(rbind, Filter(Negate(is.null), plot_df_list))
        if (is.null(plot_df) || nrow(plot_df) == 0) return(NULL)
        
        legend_labels <- sapply(parameters, function(pred) {
            row <- result_table %>% filter(Outcome == outcome_name, Predictor == pred) %>% slice(1)
            pred_display <- ifelse(pred %in% names(display_names), display_names[pred], pred)
            if (nrow(row) == 0) return(sprintf("%s: NA", pred_display))
            
            base_label <- sprintf("%s: AUC = %.3f (%s)", pred_display, row$AUC, row$CI)
            if ("p_value" %in% names(row) && !is.na(row$p_value)) {
                p_text <- if (row$p_value < 0.001) "<0.001" else sprintf("%.3f", row$p_value)
                sprintf("%s, p=%s", base_label, p_text)
            } else base_label
        })
        
        ggplot(plot_df, aes(x = fpr, y = tpr, colour = predictor)) +
            geom_line(linewidth = 1.2) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1.1, colour = "grey50") +
            scale_colour_manual(
                name = sprintf("Indices (vs %s)", ref_label),
                values = colour_palette,
                breaks = parameters,
                labels = legend_labels
            ) +
            labs(
                title = sprintf("%s (%s)", outcome_display_names[outcome_name], cohort_label),
                x = "1 – Specificity",
                y = "Sensitivity"
            ) +
            theme_minimal(base_size = 14) +
            theme(
                plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                axis.title = element_text(size = 12, face = "bold"),
                legend.position = "bottom",
                legend.text = element_text(size = 10, hjust = 0),
                legend.key.width = unit(1, "cm"),
                plot.margin = margin(10,10,10,10)
            ) +
            guides(colour = guide_legend(ncol = 1))
    }
    
    # ----------------------------- 6) RUN EVERYTHING -----------------------------
    if (!exists("data")) stop("Thiếu object `data` trong môi trường. Gán: data <- MASLD")
    
    data_cleaned   <- clean_initial_data(data)
    data_processed <- create_features(data_cleaned)
    
    data_final <- data_processed %>% filter(!is.na(masld_group))
    cat("N eligible for analysis:", nrow(data_final), "\n")
    
    params_to_analyze <- c(
        "Alb_NLR+HSI",
        "HSI",
        "FIB4",
        "Alb_NLR",
        "SII",
        "SIRI",
        "AISI",
        "NAR",
        "PNI"
    )
    
    ref_predictor <- "Alb_NLR"
    outcomes_to_analyze <- c("masld_binary", "slf_binary")
    combos <- params_to_analyze[grepl("\\+", params_to_analyze)]
    
    # 6.1 Fit combined-score logistic models + predictions + equations + coefficients (END-TO-END)
    combo_out <- fit_combined_models_and_write(
        df = data_final,
        outcomes = outcomes_to_analyze,
        combos = combos,
        min_n = 40,
        out_eq_csv = "CombinedModels_Equations.csv",
        out_coef_csv = "CombinedModels_Coefficients.csv"
    )
    data_final <- combo_out$data
    
    # 6.2 ROC analysis
    all_results <- list()
    plots <- list()
    
    for (outcome in outcomes_to_analyze) {
        cat(sprintf("\n--- ROC analysis for outcome: %s ---\n", outcome))
        
        if (outcome == "slf_binary") {
            d_use <- data_final %>% filter(masld_group == "MASLD")
            cohort_label <- "MASLD cohort"
        } else {
            d_use <- data_final
            cohort_label <- "Eligible cohort"
        }
        
        res <- analyze_roc_with_delong(d_use, outcome, params_to_analyze, ref_predictor)
        all_results[[outcome]] <- res
        print(res)
        
        plots[[outcome]] <- plot_roc_curves_delong(
            d_use, res, outcome,
            parameters = params_to_analyze,
            cohort_label = cohort_label,
            ref_label = "Alb/NLR"
        )
    }
    
    res_all <- bind_rows(all_results, .id = "OutcomeKey") %>%
        select(-OutcomeKey) %>%
        mutate(
            p_value = if_else(is.na(p_value), NA_real_, p_value),
            p_value_txt = case_when(
                is.na(p_value) ~ NA_character_,
                p_value < 0.001 ~ "<0.001",
                TRUE ~ sprintf("%.3f", p_value)
            )
        )
    
    write.csv(res_all, "ROC_DeLong_SelectedPredictors.csv", row.names = FALSE)
    cat("--- Saved: 'ROC_DeLong_SelectedPredictors.csv' ---\n")
    
    cat("\n--- COMBINE PLOTS ---\n")
    if (!is.null(plots$masld_binary) && !is.null(plots$slf_binary)) {
        final_plot <- (plots$masld_binary) / (plots$slf_binary) +
            plot_annotation(
                title = "ROC curves (single cohort)",
                subtitle = "Selected indices compared against Alb/NLR using DeLong test",
                theme = theme(
                    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 14, hjust = 0.5)
                )
            )
        
        print(final_plot)
        ggsave("ROC_SelectedPredictors_DeLong.png", plot = final_plot, width = 16, height = 14, dpi = 300)
        cat("--- Saved: 'ROC_SelectedPredictors_DeLong.png' ---\n")
    } else {
        cat("Cannot create combined plot (one or more outcome plots missing).\n")
        if (!is.null(plots$masld_binary)) print(plots$masld_binary)
        if (!is.null(plots$slf_binary)) print(plots$slf_binary)
    }
    
    cat("\n--- DONE ---\n")
    
    ```
    
- **DCA**
    
    ```r
    # ==============================================================================
    # END-TO-END (BMJ THEME): DCA ONLY (ANLR + Treat All/None)
    # Fix: bring legend back (bottom), still:
    #  - no title/tag
    #  - no border/frame
    #  - thinner/lighter lines
    #  - dashed major gridlines (x/y)
    # Output: 2 PNG files (MASLD, SLF)
    # ==============================================================================
    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(dplyr, tidyr, ggplot2, scales)
    
    if (!requireNamespace("ggsci", quietly = TRUE)) install.packages("ggsci")
    library(ggsci)
    
    if (!exists("data")) stop("Thiếu object `data`. Gán: data <- MASLD")
    
    # ----------------------------- USER CONTROLS -----------------------------
    TH_MIN <- 0.05
    TH_MAX <- 0.60
    TH_BY  <- 0.01
    
    # ----------------------------- BMJ THEME (NO FRAME, DASH GRID, LEGEND ON) -----------------------------
    theme_BMJ <- function(base_size = 16, base_family = "sans") {
      theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
          # grid: dashed x/y major lines
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey88", linewidth = 0.55, linetype = "dashed"),
    
          # remove frame/border
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
    
          # axes
          axis.title = element_text(face = "bold", color = "black"),
          axis.text  = element_text(face = "bold", color = "black"),
          axis.ticks = element_line(color = "grey40", linewidth = 0.5),
          axis.line  = element_line(color = "grey40", linewidth = 0.6),
    
          # remove any titles
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
    
          # legend back
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text  = element_text(face = "bold", color = "black"),
          legend.key = element_blank(),
    
          plot.margin = margin(10, 12, 10, 12)
        )
    }
    
    # ----------------------------- HELPERS -----------------------------
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    ensure_cols <- function(df, cols, default = NA) {
      miss <- setdiff(cols, names(df))
      if (length(miss)) for (m in miss) df[[m]] <- default
      df
    }
    
    # ----------------------------- ALCOHOL -----------------------------
    calculate_alcohol_consumption <- function(df) {
      df <- ensure_cols(df, c("ALQ121","ALQ130"), default = NA)
      df %>%
        mutate(
          ALQ121 = safe_num(ALQ121),
          ALQ130 = safe_num(ALQ130),
          days_per_week_alcohol = case_when(
            is.na(ALQ121) ~ NA_real_,
            ALQ121 == 0 ~ 0,
            ALQ121 == 1 ~ 7,
            ALQ121 == 2 ~ 6,
            ALQ121 == 3 ~ 3.5,
            ALQ121 == 4 ~ 2,
            ALQ121 == 5 ~ 1,
            ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
            ALQ121 == 7 ~ 1 / (30.4375 / 7),
            ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
            ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
            ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
            TRUE ~ NA_real_
          ),
          ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
          weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
        )
    }
    
    # ----------------------------- CLEAN / EXCLUSION -----------------------------
    clean_initial_data <- function(df) {
      cat("\n--- START CLEANING ---\n")
      n0 <- nrow(df)
    
      df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
    
      df <- ensure_cols(
        df,
        c("RIDAGEYR","RIAGENDR","LUXSMED","LUXCAPM",
          "RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV",
          "ALQ121","ALQ130",
          "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"),
        default = NA
      )
    
      cat("1) Viral hepatitis exclusion...\n")
      b <- nrow(df)
      hep_vars <- c("RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV")
      hep_vars_exist <- intersect(hep_vars, names(df))
      if (length(hep_vars_exist) > 0) {
        df <- df %>%
          mutate(across(all_of(hep_vars_exist), safe_num)) %>%
          filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
      }
      cat("-> Removed:", b - nrow(df), "| Remaining:", nrow(df), "\n\n")
    
      cat("2) Require LSM (LUXSMED) + CAP (LUXCAPM)...\n")
      b <- nrow(df)
      df <- df %>% drop_na(LUXSMED, LUXCAPM)
      cat("-> Removed:", b - nrow(df), "| Remaining:", nrow(df), "\n\n")
    
      cat("3) Require blood for ANLR (Albumin + Neut/Lymph counts OR %)...\n")
      b <- nrow(df)
      df <- df %>%
        mutate(
          LBXSAL   = safe_num(LBXSAL),
          LBDSALSI = safe_num(LBDSALSI),
          LBDNENO  = safe_num(LBDNENO),
          LBDLYMNO = safe_num(LBDLYMNO),
          LBXNEPCT = safe_num(LBXNEPCT),
          LBXLYPCT = safe_num(LBXLYPCT),
          Albumin_g_dL_tmp = case_when(
            !is.na(LBXSAL) ~ LBXSAL,
            is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
            TRUE ~ NA_real_
          ),
          has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO) & LBDLYMNO > 0,
          has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT) & LBXLYPCT > 0
        ) %>%
        filter(!is.na(Albumin_g_dL_tmp)) %>%
        filter(has_counts | has_pct)
      cat("-> Removed:", b - nrow(df), "| Remaining:", nrow(df), "\n\n")
    
      cat("4) Alcohol exclusion...\n")
      df <- calculate_alcohol_consumption(df)
      b <- nrow(df)
      df <- df %>%
        mutate(RIAGENDR = safe_num(RIAGENDR)) %>%
        filter(
          is.na(weekly_alcohol_grams) |
            (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
            (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
        )
      cat("-> Removed:", b - nrow(df), "| Remaining:", nrow(df), "\n\n")
    
      cat("5) Age >= 18...\n")
      b <- nrow(df)
      df <- df %>%
        mutate(RIDAGEYR = safe_num(RIDAGEYR)) %>%
        filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
      cat("-> Removed:", b - nrow(df), "| Remaining:", nrow(df), "\n\n")
    
      cat("--- DONE CLEANING. Total removed:", n0 - nrow(df), "| Remaining:", nrow(df), "---\n\n")
      df
    }
    
    # ----------------------------- FEATURES (ANLR) -----------------------------
    create_features <- function(df) {
      cat("--- START FEATURE ENGINEERING ---\n")
    
      df <- ensure_cols(
        df,
        c("BMXBMI","BMXWAIST",
          "LBXGLU","LBXGH","LBXSGL",
          "DIQ010","DIQ050","DIQ070",
          "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
          "BPQ040A","BPQ090D",
          "LBDTRSI","LBDHDDSI",
          "LBXSATSI","LBXSASSI","LBXPLTSI"),
        default = NA
      )
    
      out <- df %>%
        mutate(
          RIAGENDR = safe_num(RIAGENDR),
          RIDAGEYR = safe_num(RIDAGEYR),
          LUXCAPM  = safe_num(LUXCAPM),
          LUXSMED  = safe_num(LUXSMED),
    
          BMXBMI   = safe_num(BMXBMI),
          BMXWAIST = safe_num(BMXWAIST),
    
          LBXGLU   = safe_num(LBXGLU),
          LBXGH    = safe_num(LBXGH),
          LBXSGL   = safe_num(LBXSGL),
    
          DIQ010   = safe_num(DIQ010),
          DIQ050   = safe_num(DIQ050),
          DIQ070   = safe_num(DIQ070),
    
          BPXOSY1  = safe_num(BPXOSY1),
          BPXODI1  = safe_num(BPXODI1),
          BPXOSY2  = safe_num(BPXOSY2),
          BPXODI2  = safe_num(BPXODI2),
          BPXOSY3  = safe_num(BPXOSY3),
          BPXODI3  = safe_num(BPXODI3),
    
          BPQ040A  = safe_num(BPQ040A),
          BPQ090D  = safe_num(BPQ090D),
    
          LBDTRSI  = safe_num(LBDTRSI),
          LBDHDDSI = safe_num(LBDHDDSI),
    
          LBXSATSI = safe_num(LBXSATSI),
          LBXSASSI = safe_num(LBXSASSI),
          LBXPLTSI = safe_num(LBXPLTSI),
    
          LBXLYPCT = safe_num(LBXLYPCT),
          LBXNEPCT = safe_num(LBXNEPCT),
          LBDLYMNO = safe_num(LBDLYMNO),
          LBDNENO  = safe_num(LBDNENO),
          LBXSAL   = safe_num(LBXSAL),
          LBDSALSI = safe_num(LBDSALSI)
        ) %>%
        mutate(
          has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
    
          risk_bmi_waist = if_else(
            is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)), NA,
            (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)
          ),
    
          risk_glucose_diabetes = if_else(
            is.na(LBXGLU) & is.na(LBXGH) &
              (is.na(DIQ010) | DIQ010 %in% c(7,9)) &
              (is.na(DIQ050) | DIQ050 %in% c(7,9)) &
              (is.na(DIQ070) | DIQ070 %in% c(7,9)), NA,
            (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
              (if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
              (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
              (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
              (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
          ),
    
          any_bp_high =
            (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
            (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
            (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
    
          risk_blood_pressure = if_else(
            is.na(any_bp_high) & (is.na(BPQ040A) | BPQ040A %in% c(7,9)), NA,
            any_bp_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
          ),
    
          risk_triglycerides = if_else(
            is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9)), NA,
            (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          ),
    
          risk_low_hdl = if_else(
            is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9))), NA,
            ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
               (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          )
        ) %>%
        rowwise() %>%
        mutate(num_cardiometabolic_risks = sum(c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure,
                                                risk_triglycerides, risk_low_hdl), na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(
          has_one_plus_cardiometabolic_risk = if_else(
            is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
              is.na(risk_triglycerides) & is.na(risk_low_hdl), NA,
            num_cardiometabolic_risks > 0
          ),
    
          masld_group = case_when(
            has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
            has_hepatic_steatosis == FALSE ~ "non-MASLD",
            has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
            TRUE ~ NA_character_
          ),
          masld_binary = if_else(masld_group == "MASLD", 1, 0, missing = NA_real_),
    
          slf_binary = case_when(
            masld_group == "MASLD" & LUXSMED >= 8.0 ~ 1,
            masld_group == "MASLD" & LUXSMED < 8.0  ~ 0,
            TRUE ~ NA_real_
          ),
    
          Albumin_g_dL = case_when(
            !is.na(LBXSAL) ~ LBXSAL,
            is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
            TRUE ~ NA_real_
          ),
    
          Neutrophil_count   = LBDNENO,
          Lymphocyte_count   = LBDLYMNO,
          Neutrophil_percent = LBXNEPCT,
          Lymphocyte_percent = LBXLYPCT,
    
          NLR = case_when(
            !is.na(Neutrophil_count) & !is.na(Lymphocyte_count) & Lymphocyte_count > 0 ~
              Neutrophil_count / Lymphocyte_count,
            (is.na(Neutrophil_count) | is.na(Lymphocyte_count)) &
              !is.na(Neutrophil_percent) & !is.na(Lymphocyte_percent) & Lymphocyte_percent > 0 ~
              Neutrophil_percent / Lymphocyte_percent,
            TRUE ~ NA_real_
          ),
    
          ANLR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                         Albumin_g_dL / NLR, NA_real_)
        )
    
      cat("--- DONE FEATURE ENGINEERING ---\n\n")
      out
    }
    
    # ----------------------------- MODEL PROBABILITIES -----------------------------
    fit_predict_probs <- function(df, outcome_col, rhs, min_n = 80) {
      f <- as.formula(paste(outcome_col, "~", rhs))
      vars <- all.vars(f)
    
      ok <- complete.cases(df[, vars, drop = FALSE]) & df[[outcome_col]] %in% c(0,1)
      if (sum(ok) < min_n || dplyr::n_distinct(df[[outcome_col]][ok]) < 2) return(rep(NA_real_, nrow(df)))
    
      m <- tryCatch(glm(f, data = df[ok, , drop = FALSE], family = binomial()),
                    error = function(e) NULL)
      if (is.null(m)) return(rep(NA_real_, nrow(df)))
    
      p <- tryCatch(predict(m, newdata = df, type = "response"),
                    error = function(e) rep(NA_real_, nrow(df)))
      as.numeric(p)
    }
    
    # ----------------------------- DCA -----------------------------
    calculate_net_benefit <- function(y, probs, thresholds) {
      N <- length(y)
      sapply(thresholds, function(p_t) {
        tp <- sum(probs >= p_t & y == 1, na.rm = TRUE)
        fp <- sum(probs >= p_t & y == 0, na.rm = TRUE)
        nb <- (tp / N) - (fp / N) * (p_t / (1 - p_t))
        if (!is.finite(nb)) nb <- NA_real_
        nb
      })
    }
    
    generate_dca_anlr <- function(df, outcome_col, thresholds, min_n = 80) {
      probs <- fit_predict_probs(df, outcome_col, "ANLR", min_n = min_n)
      ok <- !is.na(df[[outcome_col]]) & df[[outcome_col]] %in% c(0,1) &
        !is.na(probs) & is.finite(probs)
      if (sum(ok) < min_n || dplyr::n_distinct(df[[outcome_col]][ok]) < 2) return(data.frame())
    
      y <- df[[outcome_col]][ok]
      p <- probs[ok]
      nb <- calculate_net_benefit(y, p, thresholds)
    
      data.frame(threshold = thresholds, net_benefit = nb, curve = "ANLR")
    }
    
    generate_dca_reference <- function(df, outcome_col, thresholds, min_n = 80) {
      ok <- !is.na(df[[outcome_col]]) & df[[outcome_col]] %in% c(0,1)
      if (sum(ok) < min_n || dplyr::n_distinct(df[[outcome_col]][ok]) < 2) return(data.frame())
    
      y <- df[[outcome_col]][ok]
      prev <- mean(y == 1, na.rm = TRUE)
    
      nb_all  <- prev - (1 - prev) * (thresholds / (1 - thresholds))
      nb_none <- rep(0, length(thresholds))
    
      rbind(
        data.frame(threshold = thresholds, net_benefit = nb_all,  curve = "Treat All"),
        data.frame(threshold = thresholds, net_benefit = nb_none, curve = "Treat None")
      )
    }
    
    plot_single_dca <- function(dca_anlr, dca_ref, th_min, th_max) {
      dfp <- bind_rows(dca_anlr, dca_ref) %>%
        mutate(curve = factor(curve, levels = c("ANLR", "Treat All", "Treat None")))
    
      bmj_blue <- ggsci::pal_bmj("default")(3)[1]
    
      col_map <- c("ANLR" = bmj_blue, "Treat All" = "grey70", "Treat None" = "grey70")
      lt_map  <- c("ANLR" = "solid",  "Treat All" = "dotted", "Treat None" = "solid")
      sz_map  <- c("ANLR" = 1.05,     "Treat All" = 0.80,     "Treat None" = 0.80)
    
      y_pos_max <- suppressWarnings(max(dfp$net_benefit[dfp$net_benefit >= 0], na.rm = TRUE))
      if (!is.finite(y_pos_max) || y_pos_max <= 0) y_pos_max <- 0.05
      y_top <- y_pos_max * 1.05
    
      ggplot(dfp, aes(x = threshold, y = net_benefit, color = curve, linetype = curve, size = curve)) +
        geom_hline(yintercept = 0, color = "grey75", linewidth = 0.55) +
        geom_line(na.rm = TRUE, alpha = 0.9) +
        scale_color_manual(values = col_map) +
        scale_linetype_manual(values = lt_map) +
        scale_size_manual(values = sz_map) +
        scale_x_continuous(
          limits = c(th_min, th_max),
          breaks = seq(th_min, th_max, by = 0.1),
          labels = label_number(accuracy = 0.01)
        ) +
        scale_y_continuous(labels = label_number(accuracy = 0.01)) +
        coord_cartesian(xlim = c(th_min, th_max), ylim = c(0, y_top)) +
        labs(x = "Threshold probability", y = "Net benefit") +
        theme_BMJ() +
        guides(
          color = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(linewidth = 1.3)),
          linetype = "none",
          size = "none"
        )
    }
    
    # ----------------------------- RUN -----------------------------
    thresholds <- seq(TH_MIN, TH_MAX, by = TH_BY)
    
    data_cleaned   <- clean_initial_data(data)
    data_processed <- create_features(data_cleaned)
    
    # MASLD analysis dataset
    data_analysis <- data_processed %>%
      filter(!is.na(masld_group), !is.na(masld_binary), !is.na(ANLR))
    
    # SLF analysis dataset (within MASLD only)
    data_analysis_masld <- data_analysis %>%
      filter(masld_group == "MASLD", !is.na(slf_binary))
    
    # --- 1) MASLD DCA
    dca_anlr_masld <- generate_dca_anlr(data_analysis, "masld_binary", thresholds)
    dca_ref_masld  <- generate_dca_reference(data_analysis, "masld_binary", thresholds)
    
    p_masld <- plot_single_dca(dca_anlr_masld, dca_ref_masld, TH_MIN, TH_MAX)
    print(p_masld)
    ggsave("BMJ_DCA_ANLR_MASLD.png", plot = p_masld, width = 8.2, height = 6.2, dpi = 400)
    
    # --- 2) SLF DCA
    dca_anlr_slf <- generate_dca_anlr(data_analysis_masld, "slf_binary", thresholds)
    dca_ref_slf  <- generate_dca_reference(data_analysis_masld, "slf_binary", thresholds)
    
    p_slf <- plot_single_dca(dca_anlr_slf, dca_ref_slf, TH_MIN, TH_MAX)
    print(p_slf)
    ggsave("BMJ_DCA_ANLR_SLF.png", plot = p_slf, width = 8.2, height = 6.2, dpi = 400)
    
    cat("\n--- DONE. Files saved: BMJ_DCA_ANLR_MASLD.png | BMJ_DCA_ANLR_SLF.png ---\n")
    
    ```
    
- **CIC**
    
    ```r
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(rlang)
        library(rpart)
        library(rmda)
        library(ggsci)
    })
    
    # ==============================================================================
    # CONFIGURATION & STYLE
    # ==============================================================================
    
    # BMJ Color Palette (Base R in rmda)
    bmj_cols <- c("#0072B2", "#D55E00")
    
    # ==============================================================================
    # 1) HELPERS & CLEANING
    # ==============================================================================
    
    safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    ensure_cols <- function(df, cols, default = NA) {
        miss <- setdiff(cols, names(df))
        if (length(miss)) for (m in miss) df[[m]] <- default
        df
    }
    
    calculate_alcohol_consumption <- function(data) {
        data <- ensure_cols(data, c("ALQ121","ALQ130"), default = NA)
        data %>%
            mutate(
                ALQ121 = safe_num(ALQ121),
                ALQ130 = safe_num(ALQ130),
                days_per_week_alcohol = case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    apply_exclusion_criteria <- function(data) {
        cat("\n--- START CLEANING ---\n")
        n0 <- nrow(data)
        
        # Essential columns
        data <- ensure_cols(
            data,
            c("RIAGENDR","RIDAGEYR","LUXSMED","LUXCAPM",
              "ALQ121","ALQ130",
              "RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV",
              "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"),
            default = NA
        )
        
        # Viral hepatitis exclusion (drop if any marker == 1)
        hep_vars <- c("RHD143","LBXHBC","LBDHBG","LBXHCR","LBDHCV")
        hep_vars_exist <- intersect(hep_vars, names(data))
        if (length(hep_vars_exist) > 0) {
            data <- data %>%
                mutate(across(all_of(hep_vars_exist), safe_num)) %>%
                filter(rowSums(across(all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
        }
        
        # Require LSM + CAP
        data <- data %>% drop_na(LUXSMED, LUXCAPM)
        
        # STRICT BLOOD REQUIREMENT (Albumin + Neut/Lymph counts OR %)
        data <- data %>%
            mutate(
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXLYPCT = safe_num(LBXLYPCT),
                
                Albumin_g_dL_tmp = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                has_counts = !is.na(LBDNENO) & !is.na(LBDLYMNO) & LBDLYMNO > 0,
                has_pct    = !is.na(LBXNEPCT) & !is.na(LBXLYPCT) & LBXLYPCT > 0
            ) %>%
            filter(!is.na(Albumin_g_dL_tmp)) %>%
            filter(has_counts | has_pct)
        
        # Alcohol exclusion
        data <- calculate_alcohol_consumption(data) %>%
            mutate(RIAGENDR = safe_num(RIAGENDR)) %>%
            filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        
        # Age
        data <- data %>%
            mutate(RIDAGEYR = safe_num(RIDAGEYR)) %>%
            filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
        
        cat("--- DONE CLEANING. Total removed:", n0 - nrow(data), "| Remaining:", nrow(data), "---\n")
        data
    }
    
    # ==============================================================================
    # 2) FEATURE ENGINEERING (ALNR ONLY FOR PLOTTING; OTHERS KEPT FOR MASLD/SLF DEF)
    # ==============================================================================
    
    create_features <- function(data){
        data <- ensure_cols(
            data,
            c("BMXBMI","BMXWAIST","LBXGLU","LBXGH","DIQ010",
              "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
              "BPQ040A","BPQ090D",
              "LBDTRSI","LBDHDDSI",
              "LBXSATSI","LBXSASSI","LBXPLTSI",
              "LBXTR","LBXSGTSI"),
            default = NA
        )
        
        data %>%
            mutate(
                RIAGENDR = safe_num(RIAGENDR),
                RIDAGEYR = safe_num(RIDAGEYR),
                BMXBMI   = safe_num(BMXBMI),
                BMXWAIST = safe_num(BMXWAIST),
                
                LBXSATSI = safe_num(LBXSATSI),
                LBXSASSI = safe_num(LBXSASSI),
                LBXPLTSI = safe_num(LBXPLTSI),
                
                LBDTRSI  = safe_num(LBDTRSI),
                LBDHDDSI = safe_num(LBDHDDSI),
                
                LBXGLU   = safe_num(LBXGLU),
                LBXGH    = safe_num(LBXGH),
                DIQ010   = safe_num(DIQ010),
                
                BPQ040A  = safe_num(BPQ040A),
                BPQ090D  = safe_num(BPQ090D),
                
                # --- MASLD ---
                has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, safe_num(LUXCAPM) >= 263),
                
                risk_bmi_waist = (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80),
                risk_glucose   = (LBXGLU >= 100) | (LBXGH >= 5.7) | (DIQ010 == 1),
                
                mean_SBP = rowMeans(dplyr::select(., starts_with("BPXOSY")), na.rm = TRUE),
                mean_DBP = rowMeans(dplyr::select(., starts_with("BPXODI")), na.rm = TRUE),
                risk_bp  = (mean_SBP >= 130) | (mean_DBP >= 85) | (BPQ040A == 1),
                
                risk_tg  = (LBDTRSI >= 1.70) | (BPQ090D == 1),
                risk_hdl = ((RIAGENDR == 1 & LBDHDDSI < 1.0) | (RIAGENDR == 2 & LBDHDDSI < 1.3)) | (BPQ090D == 1),
                
                metabolic_risk_count = rowSums(cbind(risk_bmi_waist, risk_glucose, risk_bp, risk_tg, risk_hdl), na.rm = TRUE),
                
                masld_group  = case_when(has_hepatic_steatosis & metabolic_risk_count >= 1 ~ "MASLD", TRUE ~ "non-MASLD"),
                masld_binary = if_else(masld_group == "MASLD", 1L, 0L),
                
                slf_group = case_when(
                    masld_group == "MASLD" & safe_num(LUXSMED) >= 8.0 ~ "SLF",
                    masld_group == "MASLD" & safe_num(LUXSMED) < 8.0  ~ "non-SLF",
                    TRUE ~ NA_character_
                ),
                slf_binary = if_else(slf_group == "SLF", 1L, 0L),
                
                # --- ALNR (Albumin / NLR) ---
                LBXSAL   = safe_num(LBXSAL),
                LBDSALSI = safe_num(LBDSALSI),
                LBDNENO  = safe_num(LBDNENO),
                LBDLYMNO = safe_num(LBDLYMNO),
                LBXNEPCT = safe_num(LBXNEPCT),
                LBXLYPCT = safe_num(LBXLYPCT),
                
                Albumin_g_dL = case_when(
                    !is.na(LBXSAL) ~ LBXSAL,
                    !is.na(LBDSALSI) ~ LBDSALSI / 10,
                    TRUE ~ NA_real_
                ),
                
                Neutrophil_val = if_else(!is.na(LBDNENO), LBDNENO, LBXNEPCT),
                Lymphocyte_val = if_else(!is.na(LBDLYMNO), LBDLYMNO, LBXLYPCT),
                
                NLR  = if_else(!is.na(Neutrophil_val) & !is.na(Lymphocyte_val) & Lymphocyte_val > 0,
                               Neutrophil_val / Lymphocyte_val, NA_real_),
                
                ALNR = if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                               Albumin_g_dL / NLR, NA_real_)
            )
    }
    
    # ==============================================================================
    # 3) CLINICAL IMPACT CURVE (BMJ STYLED) — ALNR ONLY
    #    Fix: rename plot/legend to ANLR
    # ==============================================================================
    
    fit_and_plot_bmj <- function(df, outcome_var, formula_rhs, title_text, pop_size = 1000) {
        
        temp_form <- as.formula(paste(outcome_var, "~", formula_rhs))
        vars_needed <- unique(all.vars(temp_form))
        
        model_data <- df %>%
            dplyr::select(all_of(vars_needed)) %>%
            tidyr::drop_na()
        
        if (nrow(model_data) < 50 || length(unique(model_data[[outcome_var]])) < 2) {
            plot.new()
            text(0.5, 0.5, "Insufficient Data")
            title(main = title_text, col.main = "#333333", font.main = 2)
            return(invisible(NULL))
        }
        
        model_data$y_fac <- factor(model_data[[outcome_var]], levels = c(0,1), labels = c("0","1"))
        fit_form <- as.formula(paste("y_fac ~", formula_rhs))
        
        fit <- rpart::rpart(
            fit_form,
            data = model_data,
            method = "class",
            control = rpart.control(cp = 0.001, minbucket = 20)
        )
        
        # predicted risk
        model_data$risk_hat <- predict(fit, newdata = model_data, type = "prob")[,2]
        
        # legend label control: use ANLR instead of risk_hat
        model_data$ANLR <- model_data$risk_hat
        
        model_data$y_num <- as.numeric(as.character(model_data$y_fac))
        
        dca <- rmda::decision_curve(
            y_num ~ ANLR,
            data = model_data,
            fitted.risk = TRUE,
            study.design = "cohort",
            policy = "opt-in",
            bootstraps = 50
        )
        
        rmda::plot_clinical_impact(
            dca,
            population.size = pop_size,
            col = bmj_cols,
            confidence.intervals = FALSE,
            cost.benefit.axis = FALSE,
            legend.position = "topright"
        )
        
        grid(col = "lightgray", lty = "dotted")
    }
    
    # ==============================================================================
    # 4) EXECUTION (ALNR ONLY)
    # ==============================================================================
    
    if (!exists("MASLD")) stop("Please load your dataset into 'MASLD' first.")
    
    df_clean <- apply_exclusion_criteria(MASLD)
    df_feat  <- create_features(df_clean)
    
    # Figure 1: MASLD (single panel: ANLR label on plot)
    cat("\nGenerating Figure 1: MASLD Clinical Impact Curve (ANLR)...\n")
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1), family = "sans")
    fit_and_plot_bmj(df_feat, "masld_binary", "ALNR", "ANLR")
    
    # Figure 2: SLF within MASLD (single panel: ANLR label on plot)
    cat("Generating Figure 2: SLF Clinical Impact Curve (ANLR) within MASLD...\n")
    df_slf_sub <- df_feat %>% filter(masld_group == "MASLD")
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1), family = "sans")
    fit_and_plot_bmj(df_slf_sub, "slf_binary", "ALNR", "ANLR")
    
    par(mfrow = c(1,1))
    cat("\n--- END OF ANALYSIS ---\n")
    
    ```
    
- **Histogram**
    
    ```r
    # ==============================================================================
    # END-TO-END: SQUARE BMJ-LIKE BOXPLOTS (NARROW PANEL) — ANLR = Albumin/NLR
    # Updates requested:
    #   1) Fixed y-limit: 0–10
    #   2) Add axis lines for both x and y (axis.line)
    #   3) Add more left/right padding so boxes don't hug the y-side
    #   4) Remove extra top/bottom whitespace (no y expand, no plot margin)
    #   5) Increase x-axis text size
    #   6) Remove “MASLD only” from title
    # Output:
    #   - A_MASLD_ANLR_square.png
    #   - E_SLF_ANLR_square.png
    #   - ANLR_square_2col.png
    # ==============================================================================
    suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(ggplot2)
    
      if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("Missing package 'cowplot'. Install: install.packages('cowplot')")
      }
      library(cowplot)
    })
    
    if (!exists("data")) stop("Object `data` not found in environment.")
    df0 <- as.data.frame(data)
    
    # -----------------------------
    # HELPERS
    # -----------------------------
    ensure_cols <- function(df, cols){
      miss <- setdiff(cols, names(df))
      if (length(miss)) for (m in miss) df[[m]] <- NA
      df
    }
    coerce_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    
    # -----------------------------
    # EXCLUSION (same logic, no hsCRP / no uric acid variables used)
    # -----------------------------
    calculate_alcohol_consumption <- function(df) {
      df <- ensure_cols(df, c("ALQ121","ALQ130","RIAGENDR"))
      df %>%
        mutate(
          ALQ121 = coerce_num(ALQ121),
          ALQ130 = coerce_num(ALQ130),
          RIAGENDR = coerce_num(RIAGENDR),
          days_per_week_alcohol = case_when(
            is.na(ALQ121) ~ NA_real_,
            ALQ121 == 0 ~ 0,
            ALQ121 == 1 ~ 7,
            ALQ121 == 2 ~ 6,
            ALQ121 == 3 ~ 3.5,
            ALQ121 == 4 ~ 2,
            ALQ121 == 5 ~ 1,
            ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
            ALQ121 == 7 ~ 1 / (30.4375 / 7),
            ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
            ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
            ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
            TRUE ~ NA_real_
          ),
          ALQ130_cleaned = if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, ALQ130),
          weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
        )
    }
    
    apply_exclusion_criteria <- function(df) {
    
      df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
    
      # 1) Viral hepatitis evidence
      hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR")
      hep_exist <- intersect(hep_vars, names(df))
      if (length(hep_exist) > 0) {
        for (v in hep_exist) df[[v]] <- coerce_num(df[[v]])
        df <- df %>% filter(rowSums(across(all_of(hep_exist), ~ .x == 1), na.rm = TRUE) == 0)
      }
    
      # 2) Require CAP/LSM + required blood
      req_liver  <- c("LUXSMED","LUXCAPM")
      req_counts <- c("LBDLYMNO","LBDNENO")
      req_perc   <- c("LBXLYPCT","LBXNEPCT")
    
      if (!all(req_liver %in% names(df)))  stop("Missing LUXSMED and/or LUXCAPM.")
      if (!all(req_counts %in% names(df))) stop("Missing LBDLYMNO and/or LBDNENO.")
      if (!all(req_perc %in% names(df)))   stop("Missing LBXLYPCT and/or LBXNEPCT.")
      if (!("LBXSAL" %in% names(df)) && !("LBDSALSI" %in% names(df))) {
        stop("Missing both LBXSAL and LBDSALSI (albumin).")
      }
    
      for (v in unique(c(req_liver, req_counts, req_perc, intersect(c("LBXSAL","LBDSALSI"), names(df))))) {
        df[[v]] <- coerce_num(df[[v]])
      }
    
      df <- df %>% tidyr::drop_na(all_of(req_liver), all_of(req_counts), all_of(req_perc))
      df <- df %>% filter(!(is.na(LBXSAL) & is.na(LBDSALSI)))
    
      # 3) Exclude heavy alcohol
      df <- calculate_alcohol_consumption(df)
      df <- df %>%
        filter(
          is.na(weekly_alcohol_grams) |
            (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
            (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
        )
    
      # 4) Age >= 18
      if ("RIDAGEYR" %in% names(df)) {
        df$RIDAGEYR <- coerce_num(df$RIDAGEYR)
        df <- df %>% filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
      }
    
      df
    }
    
    # -----------------------------
    # FEATURES (MASLD/SLF + ANLR)
    # -----------------------------
    create_features <- function(df){
    
      df <- ensure_cols(df, c(
        "RIAGENDR","RIDAGEYR","BMXBMI","BMXWAIST",
        "LBXGLU","LBXGH",
        "DIQ010","DIQ050","DIQ070",
        "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
        "BPQ040A","LBDTRSI","BPQ090D","LBDHDDSI"
      ))
    
      num_vars <- c(
        "LUXSMED","LUXCAPM",
        "RIAGENDR","RIDAGEYR","BMXBMI","BMXWAIST",
        "LBXGLU","LBXGH",
        "DIQ010","DIQ050","DIQ070",
        "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
        "BPQ040A","LBDTRSI","BPQ090D","LBDHDDSI",
        "LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI"
      )
      for (v in intersect(num_vars, names(df))) df[[v]] <- coerce_num(df[[v]])
    
      df %>%
        mutate(
          albumin_gdl = case_when(
            !is.na(LBXSAL) ~ LBXSAL,
            is.na(LBXSAL) & !is.na(LBDSALSI) ~ LBDSALSI / 10,
            TRUE ~ NA_real_
          ),
          neutrophil_count = LBDNENO,
          lymphocyte_count = LBDLYMNO,
          NLR = if_else(!is.na(neutrophil_count) & !is.na(lymphocyte_count) & lymphocyte_count > 0,
                        neutrophil_count / lymphocyte_count, NA_real_),
          ANLR = if_else(!is.na(albumin_gdl) & !is.na(NLR) & NLR > 0,
                        albumin_gdl / NLR, NA_real_),
    
          has_hepatic_steatosis = if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
    
          risk_bmi_waist = if_else(
            is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)), NA,
            (BMXBMI >= 25) |
              (RIAGENDR == 1 & BMXWAIST >= 94) |
              (RIAGENDR == 2 & BMXWAIST >= 80)
          ),
    
          risk_glucose_diabetes = if_else(
            is.na(LBXGLU) & is.na(LBXGH) &
              (is.na(DIQ010) | DIQ010 %in% c(7,9)) &
              (is.na(DIQ050) | DIQ050 %in% c(7,9)) &
              (is.na(DIQ070) | DIQ070 %in% c(7,9)),
            NA,
            (if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
              (if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
              (if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
              (if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
              (if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
          ),
    
          any_bp_measurement_high =
            (if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
            (if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
            (if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
    
          risk_blood_pressure = if_else(
            is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7,9)),
            NA,
            any_bp_measurement_high | (if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
          ),
    
          risk_triglycerides = if_else(
            is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9)),
            NA,
            (if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          ),
    
          risk_low_hdl = if_else(
            is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9))),
            NA,
            ((RIAGENDR == 1 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
               (RIAGENDR == 2 & if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
              (if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
          )
        ) %>%
        rowwise() %>%
        mutate(
          num_cardiometabolic_risks = sum(
            c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure, risk_triglycerides, risk_low_hdl),
            na.rm = TRUE
          )
        ) %>%
        ungroup() %>%
        mutate(
          has_one_plus_cardiometabolic_risk = if_else(
            is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
              is.na(risk_triglycerides) & is.na(risk_low_hdl),
            NA, num_cardiometabolic_risks > 0
          ),
          is_light_drinker_final = case_when(
            is.na(weekly_alcohol_grams) | is.na(RIAGENDR) ~ NA,
            (RIAGENDR == 1 & weekly_alcohol_grams < 210) | (RIAGENDR == 2 & weekly_alcohol_grams < 140) ~ TRUE,
            !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
            TRUE ~ NA
          ),
          masld_group = case_when(
            is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
            has_hepatic_steatosis == FALSE | has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
            has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE & is_light_drinker_final == FALSE ~ "non-MASLD",
            has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ "MASLD",
            TRUE ~ NA_character_
          ),
          masld_binary = if_else(masld_group == "MASLD", 1L, 0L),
          slf_group = case_when(
            masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED >= 8.0 ~ "SLF",
            masld_group == "MASLD" & !is.na(LUXSMED) & LUXSMED < 8.0  ~ "non-SLF",
            TRUE ~ NA_character_
          ),
          slf_binary = if_else(slf_group == "SLF", 1L, 0L)
        )
    }
    
    # -----------------------------
    # BUILD DATASETS
    # -----------------------------
    required_blood_exact <- c("LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO","LBXSAL","LBDSALSI")
    miss_req <- setdiff(required_blood_exact, names(df0))
    if (length(miss_req)) stop("Missing required blood columns: ", paste(miss_req, collapse = ", "))
    
    data_cleaned  <- apply_exclusion_criteria(df0)
    data_featured <- create_features(data_cleaned)
    
    df_masld <- data_featured %>%
      filter(!is.na(masld_binary), is.finite(ANLR)) %>%
      mutate(masld_binary = as.integer(masld_binary))
    
    df_slf <- data_featured %>%
      filter(masld_group == "MASLD", !is.na(slf_binary), is.finite(ANLR)) %>%
      mutate(slf_binary = as.integer(slf_binary))
    
    # -----------------------------
    # PLOT SETTINGS
    # -----------------------------
    cols_fill <- c("Control" = "#1F77B4", "Case" = "#D62728")
    yl_fixed  <- c(0, 10)
    
    theme_square <- function(base_size = 12){
      theme_classic(base_size = base_size) +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
    
          axis.title.y = element_text(face = "bold", color = "black", size = base_size + 2),
          axis.text.y  = element_text(face = "bold", color = "black", size = base_size + 1),
    
          # increase x tick label size
          axis.text.x  = element_text(face = "bold", color = "black", size = base_size + 4),
    
          # add axis lines (ox/oy), keep ticks off
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.line.y = element_line(color = "black", linewidth = 0.9),
          axis.ticks  = element_blank(),
    
          # remove outer whitespace
          plot.margin = margin(0, 0, 0, 0),
    
          aspect.ratio = 1
        )
    }
    
    make_square_box <- function(d, outcome_var, title_txt){
    
      dd <- d %>%
        mutate(Group = case_when(
          .data[[outcome_var]] == 1 ~ "Case",
          .data[[outcome_var]] == 0 ~ "Control",
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(Group), is.finite(ANLR)) %>%
        mutate(Group = factor(Group, levels = c("Control","Case")))
    
      if (nrow(dd) == 0L) return(ggplot() + theme_void() + ggtitle(title_txt))
    
      st <- dd %>%
        group_by(Group) %>%
        summarise(med = median(ANLR, na.rm = TRUE), .groups = "drop")
    
      ggplot(dd, aes(x = Group, y = ANLR)) +
        geom_point(
          position = position_jitter(width = 0.08, height = 0),
          size = 0.6, alpha = 0.14, color = "black"
        ) +
        geom_boxplot(
          aes(fill = Group),
          width = 0.38,
          outlier.shape = NA,
          linewidth = 0.9,
          alpha = 0.95
        ) +
        geom_label(
          data = st,
          aes(x = Group, y = med, label = sprintf("%.1f", med), fill = Group),
          inherit.aes = FALSE,
          label.size = 0,
          color = "white",
          fontface = "bold",
          size = 4,
          label.padding = unit(0.12, "lines")
        ) +
        scale_fill_manual(values = cols_fill, drop = FALSE) +
    
        # more x padding so boxes don't hug y-side
        scale_x_discrete(expand = expansion(mult = c(0.45, 0.30))) +
    
        # fixed ylim, no extra top/bottom space
        scale_y_continuous(
          breaks = c(0, 2.5, 5, 7.5, 10),
          expand = expansion(mult = c(0, 0))
        ) +
        coord_cartesian(ylim = yl_fixed) +
    
        labs(title = title_txt, x = NULL, y = "ANLR (g/dL per unit)") +
        theme_square(base_size = 12)
    }
    
    pA <- make_square_box(df_masld, "masld_binary", "A. MASLD ~ ANLR")
    pE <- make_square_box(df_slf,   "slf_binary",   "E. SLF ~ ANLR")
    
    ggsave("A_MASLD_ANLR_square.png", pA, width = 3.2, height = 3.2, dpi = 600)
    ggsave("E_SLF_ANLR_square.png",   pE, width = 3.2, height = 3.2, dpi = 600)
    
    final_2col <- cowplot::plot_grid(pA, pE, ncol = 2, align = "hv")
    ggsave("ANLR_square_2col.png", final_2col, width = 6.6, height = 3.3, dpi = 600)
    
    print(final_2col)
    
    ```
    
- PSM có hình
    
    ```r
    # ==============================================================================
    # END-TO-END: PSM (MatchIt) + Love plot + Word report (Alb/NLR + Albumin + Neut/Lymph)
    # - REMOVE: CURI, hs-CRP (LBXHSCRP), Uric acid (LBXSUA / LBDSUASI)
    # - USE/REPORT: Alb_by_NLR, Albumin (g/dL & g/L), NLR,
    #               Neutrophil/Lymphocyte counts + percentages
    # - REQUIRED BLOOD COLS: LBXLYPCT, LBXNEPCT, LBDLYMNO, LBDNENO, (LBXSAL and/or LBDSALSI)
    # - No train/validation split
    # Outputs:
    #   - PSM_AlbNLR_MASLD_vs_nonMASLD.png + .docx
    #   - PSM_AlbNLR_SLF_vs_nonSLF_in_MASLD.png + .docx
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(stringr)
        library(rlang)
        library(officer)
        library(flextable)
        library(MatchIt)
        library(cobalt)
        library(broom)
        library(ggplot2)
    })
    
    # -------------------------
    # 0) UTILITIES
    # -------------------------
    fix_border_issues <- function(ft){
        ft <- flextable::border_remove(ft)
        ft <- flextable::hline(ft, border = fp_border(color = "black", width = 1), part = "header")
        ft <- flextable::hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "body")
        ft
    }
    
    format_p_value <- function(p) {
        dplyr::case_when(
            is.na(p) ~ "NA",
            p < 0.001 ~ "<0.001",
            TRUE ~ format(round(p, 3), nsmall = 3)
        )
    }
    
    format_median_iqr <- function(x, digits = 2) {
        x <- x[!is.na(x)]
        if (length(x) == 0) return("—")
        q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
        paste0(format(round(q[2], digits), nsmall = digits), " (",
               format(round(q[1], digits), nsmall = digits), " - ",
               format(round(q[3], digits), nsmall = digits), ")")
    }
    
    check_required_cols <- function(df, cols, context = ""){
        missing <- setdiff(cols, names(df))
        if (length(missing)) {
            stop(
                paste0(
                    "THIẾU CỘT BẮT BUỘC", ifelse(nchar(context) > 0, paste0(" (", context, ")"), ""), ":\n- ",
                    paste(missing, collapse = "\n- ")
                ),
                call. = FALSE
            )
        }
        invisible(TRUE)
    }
    
    ensure_id <- function(df){
        if (!("ID" %in% names(df))) {
            if ("SEQN" %in% names(df)) df$ID <- df$SEQN else df$ID <- seq_len(nrow(df))
        }
        df
    }
    
    # -------------------------
    # 1) ALCOHOL (KEEP)
    # -------------------------
    calculate_alcohol_consumption <- function(data) {
        data %>%
            dplyr::mutate(
                days_per_week_alcohol = dplyr::case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = dplyr::if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, as.numeric(ALQ130)),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    # -------------------------
    # 2) BLOOD INDICES (NEW)
    # -------------------------
    derive_blood_indices <- function(df){
        df %>%
            dplyr::mutate(
                Albumin_g_dL = dplyr::case_when(
                    !is.na(LBXSAL) ~ as.numeric(LBXSAL),              # g/dL
                    !is.na(LBDSALSI) ~ as.numeric(LBDSALSI) / 10,    # g/L -> g/dL
                    TRUE ~ NA_real_
                ),
                Albumin_g_L = dplyr::case_when(
                    !is.na(LBDSALSI) ~ as.numeric(LBDSALSI),         # g/L
                    !is.na(LBXSAL) ~ as.numeric(LBXSAL) * 10,        # g/dL -> g/L
                    TRUE ~ NA_real_
                ),
                Neutrophil_count   = as.numeric(LBDNENO),
                Lymphocyte_count   = as.numeric(LBDLYMNO),
                Neutrophil_percent = as.numeric(LBXNEPCT),
                Lymphocyte_percent = as.numeric(LBXLYPCT),
                NLR = dplyr::if_else(!is.na(Neutrophil_count) & !is.na(Lymphocyte_count) & Lymphocyte_count > 0,
                                     Neutrophil_count / Lymphocyte_count, NA_real_),
                Alb_by_NLR = dplyr::if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                                            Albumin_g_dL / NLR, NA_real_)
            )
    }
    
    # -------------------------
    # 3) EXCLUSION (UPDATED: remove hsCRP/uric/CURI)
    # -------------------------
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG TIÊU CHÍ LOẠI TRỪ (Alb/NLR pipeline) ---\n")
        initial_rows_total <- nrow(data)
        
        # Required
        required_core <- c("LUXSMED", "LUXCAPM", "LBXLYPCT", "LBXNEPCT", "LBDLYMNO", "LBDNENO")
        check_required_cols(data, required_core, context = "CORE (FibroScan/CAP + neut/lymph)")
        if (!("LBXSAL" %in% names(data)) && !("LBDSALSI" %in% names(data))) {
            stop("THIẾU CỘT ALBUMIN: cần LBXSAL (g/dL) và/hoặc LBDSALSI (g/L).", call. = FALSE)
        }
        
        # 1) Viral hepatitis removal
        cat("1) Loại bỏ viêm gan virus...\n")
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR")
        hep_vars_exist <- intersect(hep_vars, names(data))
        if (length(hep_vars_exist) > 0) {
            data2 <- data %>%
                dplyr::filter(rowSums(dplyr::across(dplyr::all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
            cat("-> Loại:", nrow(data) - nrow(data2), "\n\n")
            data <- data2
        } else {
            cat("-> Không có cột virus hep: bỏ qua.\n\n")
        }
        
        # 2) Drop NA for FibroScan/CAP + neut/lymph + albumin (either)
        cat("2) Loại NA FibroScan/CAP + blood indices + Albumin...\n")
        before <- nrow(data)
        
        data <- data %>%
            tidyr::drop_na(dplyr::all_of(c("LUXSMED","LUXCAPM","LBXLYPCT","LBXNEPCT","LBDLYMNO","LBDNENO"))) %>%
            dplyr::filter(!(is.na(LBXSAL) & is.na(LBDSALSI)))
        
        cat("-> Loại:", before - nrow(data), "\n\n")
        
        # 3) Alcohol compute + heavy drinking exclude
        cat("3) Lọc rượu bia...\n")
        data <- calculate_alcohol_consumption(data)
        before <- nrow(data)
        data <- data %>%
            dplyr::filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        cat("-> Loại:", before - nrow(data), "\n\n")
        
        # 4) Age >= 18
        cat("4) Lọc tuổi >= 18...\n")
        if ("RIDAGEYR" %in% names(data)) {
            before <- nrow(data)
            data <- data %>% dplyr::filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
            cat("-> Loại:", before - nrow(data), "\n\n")
        } else {
            cat("-> Không có RIDAGEYR: bỏ qua.\n\n")
        }
        
        cat("--- HOÀN TẤT LOẠI TRỪ: loại tổng ", initial_rows_total - nrow(data),
            " | còn ", nrow(data), " ---\n\n")
        data
    }
    
    # -------------------------
    # 4) DIAGNOSTIC CRITERIA (NO CURI)
    # -------------------------
    define_diagnostic_criteria <- function(data) {
        cat("Bắt đầu: MASLD + SLF (NO CURI)\n")
        
        data_diag <- data %>%
            dplyr::mutate(has_hepatic_steatosis = dplyr::if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263)) %>%
            dplyr::mutate(
                risk_bmi_waist = dplyr::if_else(
                    is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                    NA,
                    (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)
                ),
                risk_glucose_diabetes = dplyr::if_else(
                    is.na(LBXGLU) & is.na(LBXGH) &
                        (is.na(DIQ010) | DIQ010 %in% c(7, 9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7, 9)) &
                        (is.na(DIQ070) | DIQ070 %in% c(7, 9)),
                    NA,
                    (dplyr::if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (dplyr::if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
                        (dplyr::if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (dplyr::if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (dplyr::if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                any_bp_measurement_high =
                    (dplyr::if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (dplyr::if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (dplyr::if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                risk_blood_pressure = dplyr::if_else(
                    is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7, 9)),
                    NA,
                    any_bp_measurement_high | (dplyr::if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                risk_triglycerides = dplyr::if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9)),
                    NA,
                    (dplyr::if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (dplyr::if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                risk_low_hdl = dplyr::if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7, 9))),
                    NA,
                    ((RIAGENDR == 1 & dplyr::if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & dplyr::if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (dplyr::if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                )
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                num_cardiometabolic_risks = sum(c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure,
                                                  risk_triglycerides, risk_low_hdl), na.rm = TRUE),
                has_one_plus_cardiometabolic_risk = dplyr::if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA,
                    num_cardiometabolic_risks > 0
                )
            ) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                is_light_drinker_final = dplyr::case_when(
                    is.na(weekly_alcohol_grams) ~ NA,
                    RIAGENDR == 1 & weekly_alcohol_grams < 210 ~ TRUE,
                    RIAGENDR == 2 & weekly_alcohol_grams < 140 ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                masld_group = dplyr::case_when(
                    has_hepatic_steatosis == FALSE ~ "non-MASLD",
                    has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ dplyr::case_when(
                        is.na(is_light_drinker_final) ~ "MASLD",
                        is_light_drinker_final == TRUE ~ "MASLD",
                        is_light_drinker_final == FALSE ~ "non-MASLD",
                        TRUE ~ NA_character_
                    ),
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    TRUE ~ NA_character_
                ),
                slf_group = dplyr::case_when(
                    masld_group != "MASLD" ~ NA_character_,
                    is.na(LUXSMED) ~ NA_character_,
                    LUXSMED >= 8.0 ~ "SLF",
                    LUXSMED < 8.0 ~ "non-SLF",
                    TRUE ~ NA_character_
                )
            )
        
        cat("Hoàn tất: MASLD + SLF\n\n")
        data_diag
    }
    
    # -------------------------
    # 5) COVARIATES (PSM)
    # -------------------------
    define_categorical_variables <- function(data) {
        cat("Bắt đầu: covariates\n")
        
        data_cat <- data %>%
            dplyr::mutate(
                Age = as.numeric(RIDAGEYR),
                
                creatinine_mg_dl = dplyr::case_when(
                    "LBXSCR" %in% names(.) ~ as.numeric(LBXSCR),
                    "LBDSCRSI" %in% names(.) ~ as.numeric(LBDSCRSI) / 88.4,
                    TRUE ~ NA_real_
                ),
                eGFR = dplyr::if_else(
                    !is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                    175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                        dplyr::if_else(RIAGENDR == 2, 0.742, 1) * dplyr::if_else(RIDRETH1 == 4, 1.212, 1),
                    NA_real_
                ),
                
                Gender = factor(dplyr::case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female"),
                                levels = c("Male", "Female")),
                
                Education_Level = factor(dplyr::case_when(
                    DMDEDUC2 == 1 ~ "Less than 9th grade",
                    DMDEDUC2 == 2 ~ "9-11th grade",
                    DMDEDUC2 == 3 ~ "High school graduate/GED",
                    DMDEDUC2 == 4 ~ "Some college or AA degree",
                    DMDEDUC2 == 5 ~ "College graduate or above",
                    TRUE ~ NA_character_
                ), levels = c("Less than 9th grade","9-11th grade","High school graduate/GED",
                              "Some college or AA degree","College graduate or above")),
                
                Marital_Status = factor(dplyr::case_when(
                    DMDMARTZ == 1 ~ "Married",
                    DMDMARTZ == 2 ~ "Widowed",
                    DMDMARTZ == 3 ~ "Divorced",
                    DMDMARTZ == 4 ~ "Separated",
                    DMDMARTZ == 5 ~ "Never married",
                    DMDMARTZ == 6 ~ "Living with partner",
                    TRUE ~ NA_character_
                ), levels = c("Married","Living with partner","Never married","Divorced","Widowed","Separated")),
                
                Race_Ethnicity = factor(dplyr::case_when(
                    RIDRETH1 == 1 ~ "Mexican American",
                    RIDRETH1 == 2 ~ "Other Hispanic",
                    RIDRETH1 == 3 ~ "Non-Hispanic White",
                    RIDRETH1 == 4 ~ "Non-Hispanic Black",
                    RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
                    TRUE ~ NA_character_
                )),
                
                Smoking_Status = factor(dplyr::case_when(
                    SMQ020 == 2 ~ "Non-smoker",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
                    TRUE ~ NA_character_
                ), levels = c("Non-smoker","Former smoker","Current smoker")),
                
                Obesity = factor(dplyr::if_else(!is.na(BMXBMI) & BMXBMI >= 30, "Yes", "No", missing = "No"),
                                 levels = c("No","Yes")),
                
                Diabetes_Mellitus = factor(dplyr::if_else(
                    DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                        (!is.na(LBXGH) & LBXGH >= 6.5) |
                        (!is.na(LBXSGL) & LBXSGL >= 126),
                    "Yes", "No", missing = "No"
                ), levels = c("No","Yes")),
                
                Hypertension = factor(dplyr::if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"),
                                      levels = c("No","Yes")),
                
                Chronic_Kidney_Disease = factor(dplyr::if_else(
                    (KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60)),
                    "Yes", "No", missing = "No"
                ), levels = c("No","Yes")),
                
                Cancer = factor(dplyr::if_else(MCQ220 == 1, "Yes", "No", missing = "No"),
                                levels = c("No","Yes")),
                
                Cardiovascular_Disease = factor(dplyr::if_else(
                    MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1,
                    "Yes", "No", missing = "No"
                ), levels = c("No","Yes"))
            )
        
        cat("Hoàn tất: covariates\n\n")
        data_cat
    }
    
    # -------------------------
    # 6) FEATURE PIPELINE
    # -------------------------
    create_all_features <- function(data) {
        cat("--- BẮT ĐẦU: FEATURE ENGINEERING ---\n")
        
        optional_cols <- c(
            "LBXSGL","SMQ040","DIQ050","DIQ070","BPQ040A","KIQ022","KIQ025",
            "MCQ220","MCQ160B","MCQ160C","MCQ160D","MCQ160E","BPQ020",
            "BMXBMI","BMXWAIST","LBXGLU","LBXGH","DMDEDUC2","DMDMARTZ",
            "RIDRETH1","SMQ020","BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3",
            "BPQ090D","LBDTRSI","LBDHDDSI","ALQ121","ALQ130"
        )
        for (col in optional_cols) if (!col %in% names(data)) data[[col]] <- NA
        
        data %>%
            ensure_id() %>%
            derive_blood_indices() %>%
            define_diagnostic_criteria() %>%
            define_categorical_variables() %>%
            { cat("--- HOÀN TẤT: FEATURE ENGINEERING ---\n\n"); . }
    }
    
    # -------------------------
    # 7) BALANCE TABLE (FIXED namespaces)
    # -------------------------
    create_balance_table <- function(psm_model) {
        cat("Tạo bảng cân bằng covariates...\n")
        
        bal_data <- cobalt::bal.tab(psm_model, binary = "std", abs = TRUE, un = TRUE)
        balance_df <- as.data.frame(bal_data$Balance)
        balance_df$Variable <- rownames(balance_df)
        
        balance_df_formatted <- balance_df %>%
            dplyr::select(Variable, Diff.Unmatched, Diff.Matched) %>%
            dplyr::rename(
                `Biến số` = Variable,
                `SMD chưa ghép` = Diff.Unmatched,
                `SMD đã ghép` = Diff.Matched
            ) %>%
            dplyr::mutate(dplyr::across(where(is.numeric), ~ format(round(., 3), nsmall = 3)))
        
        flextable::flextable(balance_df_formatted) %>%
            flextable::autofit() %>%
            flextable::theme_booktabs() %>%
            flextable::set_caption("Bảng 1: Cân bằng biến số trước và sau ghép (SMD)") %>%
            flextable::align(j = 2:3, align = "center", part = "all") %>%
            fix_border_issues()
    }
    
    # -------------------------
    # 8) OUTCOME TABLE (ONE OUTCOME) - namespaces explicit
    # -------------------------
    create_outcome_table_one <- function(data_pre, data_post, psm_vars,
                                         outcome_var,
                                         control_level, treatment_level, Ns,
                                         analysis_caption) {
        
        unmatched_g0 <- data_pre %>% dplyr::filter(group == 0) %>% dplyr::pull(!!rlang::sym(outcome_var))
        unmatched_g1 <- data_pre %>% dplyr::filter(group == 1) %>% dplyr::pull(!!rlang::sym(outcome_var))
        
        val_g0_un <- format_median_iqr(unmatched_g0)
        val_g1_un <- format_median_iqr(unmatched_g1)
        p_un <- tryCatch(stats::wilcox.test(unmatched_g0, unmatched_g1)$p.value, error = function(e) NA_real_)
        
        dp <- data_post %>%
            dplyr::filter(!is.na(subclass)) %>%
            dplyr::group_by(subclass) %>%
            dplyr::filter(dplyr::n_distinct(group) == 2) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(subclass, group)
        
        matched_g0 <- dp %>% dplyr::filter(group == 0) %>% dplyr::pull(!!rlang::sym(outcome_var))
        matched_g1 <- dp %>% dplyr::filter(group == 1) %>% dplyr::pull(!!rlang::sym(outcome_var))
        
        val_g0_m <- format_median_iqr(matched_g0)
        val_g1_m <- format_median_iqr(matched_g1)
        p_m <- tryCatch(stats::wilcox.test(matched_g0, matched_g1, paired = TRUE)$p.value, error = function(e) NA_real_)
        
        formula_adjusted <- stats::as.formula(paste(outcome_var, "~ group +", paste(psm_vars, collapse = " + ")))
        model_adj <- stats::lm(formula_adjusted, data = data_post, weights = data_post$weights)
        
        ci_adj <- tryCatch(stats::confint(model_adj), error = function(e) NULL)
        tidy_adj <- broom::tidy(model_adj) %>% dplyr::filter(term == "group")
        
        beta_formatted <- "—"
        p_adj <- NA_real_
        if (nrow(tidy_adj) == 1 && !is.null(ci_adj) && ("group" %in% rownames(ci_adj))) {
            beta_formatted <- paste0(
                format(round(tidy_adj$estimate, 2), nsmall = 2), " (",
                format(round(ci_adj["group", 1], 2), nsmall = 2), " - ",
                format(round(ci_adj["group", 2], 2), nsmall = 2), ")"
            )
            p_adj <- tidy_adj$p.value
        }
        
        out_df <- data.frame(
            Categories = c(paste0(outcome_var, ", Median (IQR)"), "Adjusted Beta (95% CI)"),
            Original_Control = c(val_g0_un, "—"),
            Original_Treatment = c(val_g1_un, "—"),
            Original_P = c(format_p_value(p_un), "—"),
            Matched_Control = c(val_g0_m, "—"),
            Matched_Treatment = c(val_g1_m, beta_formatted),
            Matched_P = c(format_p_value(p_m), format_p_value(p_adj)),
            stringsAsFactors = FALSE
        )
        
        flextable::flextable(out_df) %>%
            flextable::set_caption(analysis_caption) %>%
            flextable::set_header_labels(
                Categories = "Categories",
                Original_Control = paste0(control_level, " (N=", Ns$n_c_un, ")"),
                Original_Treatment = paste0(treatment_level, " (N=", Ns$n_t_un, ")"),
                Original_P = "P-value",
                Matched_Control = paste0(control_level, " (N=", Ns$n_c_m, ")"),
                Matched_Treatment = paste0(treatment_level, " (N=", Ns$n_t_m, ")"),
                Matched_P = "P-value"
            ) %>%
            flextable::add_header_row(values = c("", "Original Cohort", "Matched Cohort"),
                                      colwidths = c(1, 3, 3)) %>%
            flextable::theme_booktabs() %>%
            flextable::autofit() %>%
            flextable::align(j = 2:7, align = "center", part = "all") %>%
            flextable::align(j = 1, align = "left", part = "all") %>%
            flextable::bold(i = 1:2, part = "header", bold = TRUE) %>%
            fix_border_issues() %>%
            flextable::merge_at(i = 2, j = 2:3, part = "body") %>%
            flextable::merge_at(i = 2, j = 5, part = "body") %>%
            flextable::align(i = 2, j = c(2,5,6,7), align = "center", part = "body")
    }
    
    # -------------------------
    # 9) PSM MULTI-OUTCOMES (namespaces explicit)
    # -------------------------
    run_psm_analysis_multi_outcomes <- function(data_input,
                                                group_var_string,
                                                treatment_level,
                                                control_level,
                                                psm_vars,
                                                outcomes,
                                                output_filename,
                                                table_caption_prefix = "Bảng: So sánh") {
        
        cat(paste0("\n--- BẮT ĐẦU PSM: ", group_var_string, " | ", treatment_level, " vs ", control_level, " ---\n"))
        
        data_pre <- data_input %>%
            dplyr::filter(!!rlang::sym(group_var_string) %in% c(treatment_level, control_level)) %>%
            dplyr::mutate(group = ifelse(!!rlang::sym(group_var_string) == treatment_level, 1, 0)) %>%
            dplyr::select(ID, group, dplyr::all_of(psm_vars), dplyr::all_of(outcomes), dplyr::everything()) %>%
            tidyr::drop_na(dplyr::all_of(psm_vars)) %>%
            dplyr::filter(!is.na(Alb_by_NLR))
        
        Ns0 <- sum(data_pre$group == 0)
        Ns1 <- sum(data_pre$group == 1)
        
        cat("N trước ghép:", nrow(data_pre), " | control:", Ns0, " | treatment:", Ns1, "\n")
        if (Ns0 < 10 || Ns1 < 10) {
            cat("CẢNH BÁO: Không đủ quan sát cho PSM. Bỏ qua.\n")
            return(invisible(NULL))
        }
        
        formula_psm <- stats::as.formula(paste("group ~", paste(psm_vars, collapse = " + ")))
        psm_model <- MatchIt::matchit(
            formula = formula_psm,
            data = data_pre,
            method = "nearest",
            distance = "glm",
            caliper = 0.1,
            ratio = 1,
            replace = FALSE
        )
        
        data_post <- MatchIt::match.data(psm_model, data = data_pre)
        
        Ns_list <- list(
            n_c_un = sum(data_pre$group == 0),
            n_t_un = sum(data_pre$group == 1),
            n_c_m  = sum(data_post$group == 0),
            n_t_m  = sum(data_post$group == 1)
        )
        cat("N sau ghép:", nrow(data_post), " | control:", Ns_list$n_c_m, " | treatment:", Ns_list$n_t_m, "\n")
        
        # Love plot
        cat("Tạo Love plot...\n")
        plot_filename <- gsub("\\.docx$", ".png", output_filename)
        
        variable_labels <- c(
            "Age" = "Age, years",
            "GenderMale" = "Male",
            "GenderFemale" = "Female",
            "Race_Ethnicity" = "Race/Ethnicity",
            "Education_Level" = "Education",
            "Marital_Status" = "Marital status",
            "Smoking_Status" = "Smoking",
            "ObesityYes" = "Obesity (BMI ≥ 30)",
            "Diabetes_MellitusYes" = "Diabetes mellitus",
            "HypertensionYes" = "Hypertension",
            "Chronic_Kidney_DiseaseYes" = "CKD",
            "CancerYes" = "Cancer",
            "Cardiovascular_DiseaseYes" = "CVD"
        )
        
        lp <- tryCatch({
            cobalt::love.plot(
                psm_model,
                stats = "m",
                binary = "std",
                abs = TRUE,
                thresholds = c(m = 0.2),
                var.names = variable_labels
            ) +
                ggplot2::theme_minimal(base_size = 12) +
                ggplot2::labs(
                    title = paste("Covariate Balance:", treatment_level, "vs.", control_level),
                    x = "Absolute Standardized Mean Difference",
                    y = NULL
                )
        }, error = function(e) {
            cat("-> LỖI Love plot:", e$message, "\n")
            NULL
        })
        
        if (!is.null(lp)) {
            ggplot2::ggsave(plot_filename, plot = lp, width = 8, height = 6, dpi = 300, bg = "white")
            cat("-> Đã lưu Love plot:", normalizePath(plot_filename), "\n")
        }
        
        # Word
        cat("Xuất Word...\n")
        doc <- officer::read_docx()
        
        ft_bal <- create_balance_table(psm_model)
        doc <- officer::body_add_flextable(doc, value = ft_bal)
        doc <- officer::body_add_par(doc, " ", style = "Normal")
        
        for (outcome_var in outcomes) {
            if (!outcome_var %in% names(data_pre) || !outcome_var %in% names(data_post)) next
            
            cap <- paste0(table_caption_prefix, " ", outcome_var, " giữa nhóm ",
                          treatment_level, " và ", control_level)
            
            ft_out <- create_outcome_table_one(
                data_pre  = data_pre  %>% dplyr::filter(!is.na(.data[[outcome_var]])),
                data_post = data_post %>% dplyr::filter(!is.na(.data[[outcome_var]])),
                psm_vars = psm_vars,
                outcome_var = outcome_var,
                control_level = control_level,
                treatment_level = treatment_level,
                Ns = Ns_list,
                analysis_caption = cap
            )
            doc <- officer::body_add_flextable(doc, value = ft_out)
            doc <- officer::body_add_par(doc, " ", style = "Normal")
        }
        
        note_text <- paste("Ghi chú: Mô hình 'Đã ghép, đã điều chỉnh' điều chỉnh cho:",
                           paste(psm_vars, collapse = ", "))
        doc <- officer::body_add_par(doc, note_text, style = "Normal")
        doc <- officer::body_add_par(doc, "Adjusted Beta (95% CI): khác biệt của nhóm treatment (group=1) so với control (group=0).",
                                     style = "Normal")
        
        print(doc, target = output_filename)
        cat("-> Đã lưu Word:", normalizePath(output_filename), "\n")
        
        invisible(list(psm_model = psm_model, data_pre = data_pre, data_post = data_post))
    }
    
    # ==============================================================================
    # 10) MAIN RUN (NO SPLIT)
    # ==============================================================================
    if (!exists("data") || !is.data.frame(data)) {
        stop("LỖI: Biến 'data' không tồn tại hoặc không phải data.frame.", call. = FALSE)
    }
    
    # Exclusion
    data_cleaned <- apply_exclusion_criteria(data)
    
    # Features
    data_featured <- create_all_features(data_cleaned)
    
    # Covariates for PSM (use Age created above, not RIDAGEYR)
    psm_vars_shared <- c(
        "Age",
        "Gender",
        "Race_Ethnicity",
        "Marital_Status",
        "Education_Level",
        "Smoking_Status",
        "Obesity",
        "Diabetes_Mellitus",
        "Hypertension",
        "Chronic_Kidney_Disease",
        "Cancer",
        "Cardiovascular_Disease"
    )
    
    check_required_cols(data_featured, psm_vars_shared, context = "PSM covariates")
    
    # Outcomes to report
    outcomes_report <- c(
        "Alb_by_NLR",
        "Albumin_g_dL",
        "Albumin_g_L",
        "NLR",
        "Neutrophil_count",
        "Lymphocyte_count",
        "Neutrophil_percent",
        "Lymphocyte_percent"
    )
    check_required_cols(data_featured, outcomes_report, context = "Alb/NLR outcomes")
    
    cat("\n--- BẮT ĐẦU QUY TRÌNH PSM (Alb/NLR) ---\n")
    
    # Analysis 1: MASLD vs non-MASLD
    run_psm_analysis_multi_outcomes(
        data_input = data_featured,
        group_var_string = "masld_group",
        treatment_level = "MASLD",
        control_level = "non-MASLD",
        psm_vars = psm_vars_shared,
        outcomes = outcomes_report,
        output_filename = "PSM_AlbNLR_MASLD_vs_nonMASLD.docx",
        table_caption_prefix = "Bảng: So sánh"
    )
    
    # Analysis 2: SLF vs non-SLF within MASLD
    run_psm_analysis_multi_outcomes(
        data_input = data_featured %>% dplyr::filter(masld_group == "MASLD"),
        group_var_string = "slf_group",
        treatment_level = "SLF",
        control_level = "non-SLF",
        psm_vars = psm_vars_shared,
        outcomes = outcomes_report,
        output_filename = "PSM_AlbNLR_SLF_vs_nonSLF_in_MASLD.docx",
        table_caption_prefix = "Bảng: So sánh"
    )
    
    cat("\n--- TOÀN BỘ QUY TRÌNH PSM (Alb/NLR) ĐÃ HOÀN TẤT ---\n")
    
    ```
    
- **Forrest plot + PSM**
    
    ```r
    # ==============================================================================
    # END-TO-END: PSM + Forest (3 models) — ANLR (NO CURI/hsCRP/Uric)
    # STYLE: BMJ (British Medical Journal) for Forest Plots
    # ==============================================================================
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(officer)
        library(flextable)
        library(rlang)
        library(MatchIt)
        library(cobalt)
        library(broom)
        library(ggplot2)
        library(survival)
        library(sandwich)
        library(lmtest)
        library(scales)
        library(gridExtra)
    })
    
    # -------------------------
    # 0) UTILITIES
    # -------------------------
    check_required_cols <- function(df, cols, context = "") {
        miss <- setdiff(cols, names(df))
        if (length(miss)) {
            stop(
                paste0("THIẾU CỘT BẮT BUỘC", ifelse(nchar(context) > 0, paste0(" (", context, ")"), ""), ":\n- ",
                       paste(miss, collapse = "\n- ")),
                call. = FALSE
            )
        }
        invisible(TRUE)
    }
    
    ensure_id <- function(df) {
        if (!("ID" %in% names(df))) {
            if ("SEQN" %in% names(df)) df$ID <- df$SEQN else df$ID <- seq_len(nrow(df))
        }
        df
    }
    
    safe_rename <- function(df, rename_map_new_eq_old) {
        for (new_nm in names(rename_map_new_eq_old)) {
            old_nm <- rename_map_new_eq_old[[new_nm]]
            if (old_nm %in% names(df) && !(new_nm %in% names(df))) {
                names(df)[names(df) == old_nm] <- new_nm
            }
        }
        df
    }
    
    fix_border_issues <- function(ft){
        ft <- flextable::border_remove(ft)
        ft <- flextable::hline(ft, border = fp_border(color = "black", width = 1), part = "header")
        ft <- flextable::hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "body")
        ft
    }
    
    format_p_value <- function(p) {
        dplyr::case_when(
            is.na(p) ~ "NA",
            p < 0.001 ~ "<0.001",
            TRUE ~ format(round(p, 3), nsmall = 3)
        )
    }
    
    format_median_iqr <- function(x, digits = 2) {
        x <- x[!is.na(x)]
        if (length(x) == 0) return("—")
        q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
        paste0(format(round(q[2], digits), nsmall = digits), " (",
               format(round(q[1], digits), nsmall = digits), " - ",
               format(round(q[3], digits), nsmall = digits), ")")
    }
    
    # -------------------------
    # 1) ALCOHOL & DATA PREP
    # -------------------------
    calculate_alcohol_consumption <- function(data) {
        data %>%
            dplyr::mutate(
                days_per_week_alcohol = dplyr::case_when(
                    is.na(ALQ121) ~ NA_real_,
                    ALQ121 == 0 ~ 0,
                    ALQ121 == 1 ~ 7,
                    ALQ121 == 2 ~ 6,
                    ALQ121 == 3 ~ 3.5,
                    ALQ121 == 4 ~ 2,
                    ALQ121 == 5 ~ 1,
                    ALQ121 == 6 ~ mean(c(2, 3)) / (30.4375 / 7),
                    ALQ121 == 7 ~ 1 / (30.4375 / 7),
                    ALQ121 == 8 ~ mean(c(7, 11)) / 52.1775,
                    ALQ121 == 9 ~ mean(c(3, 6)) / 52.1775,
                    ALQ121 == 10 ~ mean(c(1, 2)) / 52.1775,
                    TRUE ~ NA_real_
                ),
                ALQ130_cleaned = dplyr::if_else(ALQ130 %in% c(777, 999) | is.na(ALQ130), NA_real_, as.numeric(ALQ130)),
                weekly_alcohol_grams = ALQ130_cleaned * days_per_week_alcohol * 14
            )
    }
    
    derive_blood_indices <- function(df) {
        df %>%
            dplyr::mutate(
                Albumin_g_dL = dplyr::case_when(
                    !is.na(LBXSAL) ~ as.numeric(LBXSAL),
                    !is.na(LBDSALSI) ~ as.numeric(LBDSALSI) / 10,
                    TRUE ~ NA_real_
                ),
                Albumin_g_L = dplyr::case_when(
                    !is.na(LBDSALSI) ~ as.numeric(LBDSALSI),
                    !is.na(LBXSAL) ~ as.numeric(LBXSAL) * 10,
                    TRUE ~ NA_real_
                ),
                Neutrophil_count   = as.numeric(LBDNENO),
                Lymphocyte_count   = as.numeric(LBDLYMNO),
                Neutrophil_percent = as.numeric(LBXNEPCT),
                Lymphocyte_percent = as.numeric(LBXLYPCT),
    
                NLR = dplyr::if_else(!is.na(Neutrophil_count) & !is.na(Lymphocyte_count) & Lymphocyte_count > 0,
                                     Neutrophil_count / Lymphocyte_count, NA_real_),
    
                # ANLR = Albumin / NLR  (renamed from Alb_by_NLR)
                ANLR = dplyr::if_else(!is.na(Albumin_g_dL) & !is.na(NLR) & NLR > 0,
                                      Albumin_g_dL / NLR, NA_real_)
            )
    }
    
    # -------------------------
    # 2) EXCLUSION & DIAGNOSTICS
    # -------------------------
    apply_exclusion_criteria <- function(data) {
        cat("--- BẮT ĐẦU: ÁP DỤNG CÁC TIÊU CHÍ LOẠI TRỪ (ANLR) ---\n")
        initial_rows_total <- nrow(data)
    
        # 1) Viral hepatitis
        cat("1) Loại bỏ viêm gan virus...\n")
        hep_vars <- c("RHD143", "LBXHBC", "LBDHBG", "LBXHCR")
        hep_vars_exist <- intersect(hep_vars, names(data))
        if (length(hep_vars_exist) > 0) {
            data2 <- data %>%
                dplyr::filter(rowSums(dplyr::across(dplyr::all_of(hep_vars_exist), ~ .x == 1), na.rm = TRUE) == 0)
            cat("-> Đã loại bỏ", nrow(data) - nrow(data2), "người.\n\n")
            data <- data2
        }
    
        # 2) Required core measurements
        cat("2) Loại bỏ thiếu dữ liệu ở FibroScan/CAP + blood...\n")
        required_core <- c("LUXSMED", "LUXCAPM", "LBXLYPCT", "LBXNEPCT", "LBDLYMNO", "LBDNENO")
        check_required_cols(data, required_core, context = "CORE required")
    
        if (!("LBXSAL" %in% names(data)) && !("LBDSALSI" %in% names(data))) {
            stop("THIẾU CỘT ALBUMIN: cần LBXSAL (g/dL) hoặc LBDSALSI (g/L).", call. = FALSE)
        }
    
        before <- nrow(data)
        data <- data %>%
            tidyr::drop_na(dplyr::all_of(required_core)) %>%
            dplyr::filter(!(is.na(LBXSAL) & is.na(LBDSALSI)))
        cat("-> Đã loại bỏ", before - nrow(data), "hàng do thiếu FibroScan/CAP/blood/albumin.\n\n")
    
        # 3) Alcohol
        cat("3) Loại bỏ tiêu thụ rượu bia vượt ngưỡng...\n")
        data <- calculate_alcohol_consumption(data)
        before <- nrow(data)
        data <- data %>%
            dplyr::filter(
                is.na(weekly_alcohol_grams) |
                    (RIAGENDR == 1 & weekly_alcohol_grams <= 210) |
                    (RIAGENDR == 2 & weekly_alcohol_grams <= 140)
            )
        cat("-> Đã loại bỏ", before - nrow(data), "hàng do rượu bia.\n\n")
    
        # 4) Age
        cat("4) Loại bỏ < 18 tuổi...\n")
        if ("RIDAGEYR" %in% names(data)) {
            before <- nrow(data)
            data <- data %>% dplyr::filter(is.na(RIDAGEYR) | RIDAGEYR >= 18)
            cat("-> Đã loại bỏ", before - nrow(data), "hàng do tuổi < 18.\n\n")
        }
    
        cat("--- HOÀN TẤT LOẠI TRỪ. Còn lại:", nrow(data), "---\n\n")
        data
    }
    
    define_diagnostic_criteria <- function(data) {
        cat("Bắt đầu: Xác định MASLD + SLF (NO CURI)\n")
    
        data %>%
            dplyr::mutate(
                glucose_mg_dl = dplyr::if_else(!is.na(LBXGLU), LBXGLU, LBDGLUSI / 0.05551),
                has_hepatic_steatosis = dplyr::if_else(is.na(LUXCAPM), NA, LUXCAPM >= 263),
                risk_bmi_waist = dplyr::if_else(is.na(RIAGENDR) | (is.na(BMXBMI) & is.na(BMXWAIST)),
                                                NA,
                                                (BMXBMI >= 25) | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80)),
                risk_glucose_diabetes = dplyr::if_else(
                    is.na(LBXGLU) & is.na(LBXGH) & (is.na(DIQ010) | DIQ010 %in% c(7,9)) &
                        (is.na(DIQ050) | DIQ050 %in% c(7,9)) & (is.na(DIQ070) | DIQ070 %in% c(7,9)),
                    NA,
                    (dplyr::if_else(is.na(LBXGLU), FALSE, LBXGLU >= 100)) |
                        (dplyr::if_else(is.na(LBXGH), FALSE, LBXGH >= 5.7)) |
                        (dplyr::if_else(is.na(DIQ010), FALSE, DIQ010 == 1)) |
                        (dplyr::if_else(is.na(DIQ050), FALSE, DIQ050 == 1)) |
                        (dplyr::if_else(is.na(DIQ070), FALSE, DIQ070 == 1))
                ),
                any_bp_measurement_high =
                    (dplyr::if_else(is.na(BPXOSY1) | is.na(BPXODI1), NA, BPXOSY1 >= 130 | BPXODI1 >= 85)) |
                    (dplyr::if_else(is.na(BPXOSY2) | is.na(BPXODI2), NA, BPXOSY2 >= 130 | BPXODI2 >= 85)) |
                    (dplyr::if_else(is.na(BPXOSY3) | is.na(BPXODI3), NA, BPXOSY3 >= 130 | BPXODI3 >= 85)),
                risk_blood_pressure = dplyr::if_else(
                    is.na(any_bp_measurement_high) & (is.na(BPQ040A) | BPQ040A %in% c(7,9)),
                    NA,
                    any_bp_measurement_high | (dplyr::if_else(is.na(BPQ040A), FALSE, BPQ040A == 1))
                ),
                risk_triglycerides = dplyr::if_else(
                    is.na(LBDTRSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9)),
                    NA,
                    (dplyr::if_else(is.na(LBDTRSI), FALSE, LBDTRSI >= 1.70)) |
                        (dplyr::if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                ),
                risk_low_hdl = dplyr::if_else(
                    is.na(RIAGENDR) | (is.na(LBDHDDSI) & (is.na(BPQ090D) | BPQ090D %in% c(7,9))),
                    NA,
                    ((RIAGENDR == 1 & dplyr::if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.0)) |
                         (RIAGENDR == 2 & dplyr::if_else(is.na(LBDHDDSI), FALSE, LBDHDDSI < 1.3))) |
                        (dplyr::if_else(is.na(BPQ090D), FALSE, BPQ090D == 1))
                )
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                num_cardiometabolic_risks = sum(c(risk_bmi_waist, risk_glucose_diabetes, risk_blood_pressure,
                                                  risk_triglycerides, risk_low_hdl), na.rm = TRUE),
                has_one_plus_cardiometabolic_risk = dplyr::if_else(
                    is.na(risk_bmi_waist) & is.na(risk_glucose_diabetes) & is.na(risk_blood_pressure) &
                        is.na(risk_triglycerides) & is.na(risk_low_hdl),
                    NA,
                    num_cardiometabolic_risks > 0
                )
            ) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                is_light_drinker_final = dplyr::case_when(
                    is.na(weekly_alcohol_grams) ~ NA,
                    RIAGENDR == 1 & weekly_alcohol_grams < 210 ~ TRUE,
                    RIAGENDR == 2 & weekly_alcohol_grams < 140 ~ TRUE,
                    !is.na(weekly_alcohol_grams) & !is.na(RIAGENDR) ~ FALSE,
                    TRUE ~ NA
                ),
                masld_group = dplyr::case_when(
                    has_hepatic_steatosis == FALSE ~ "non-MASLD",
                    has_one_plus_cardiometabolic_risk == FALSE ~ "non-MASLD",
                    has_hepatic_steatosis == TRUE & has_one_plus_cardiometabolic_risk == TRUE ~ dplyr::case_when(
                        is.na(is_light_drinker_final) ~ "MASLD",
                        is_light_drinker_final == TRUE ~ "MASLD",
                        is_light_drinker_final == FALSE ~ "non-MASLD",
                        TRUE ~ NA_character_
                    ),
                    is.na(has_hepatic_steatosis) | is.na(has_one_plus_cardiometabolic_risk) ~ NA_character_,
                    TRUE ~ NA_character_
                ),
                slf_group = dplyr::case_when(
                    masld_group != "MASLD" ~ NA_character_,
                    is.na(LUXSMED) ~ NA_character_,
                    LUXSMED >= 8.0 ~ "SLF",
                    LUXSMED < 8.0 ~ "non-SLF",
                    TRUE ~ NA_character_
                )
            )
    }
    
    define_categorical_variables <- function(data) {
        cat("Bắt đầu: Tạo covariates\n")
    
        data %>%
            dplyr::mutate(
                creatinine_mg_dl = dplyr::case_when(
                    "LBXSCR" %in% names(.) ~ as.numeric(LBXSCR),
                    "LBDSCRSI" %in% names(.) ~ as.numeric(LBDSCRSI) / 88.4,
                    TRUE ~ NA_real_
                ),
                eGFR = dplyr::if_else(
                    !is.na(creatinine_mg_dl) & !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH1),
                    175 * (creatinine_mg_dl ^ -1.154) * (RIDAGEYR ^ -0.203) *
                        dplyr::if_else(RIAGENDR == 2, 0.742, 1) * dplyr::if_else(RIDRETH1 == 4, 1.212, 1),
                    NA_real_
                ),
                Age = as.numeric(RIDAGEYR),
    
                Gender = factor(dplyr::case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female"),
                                levels = c("Male", "Female")),
                Education_Level = factor(dplyr::case_when(
                    DMDEDUC2 == 1 ~ "Less than 9th grade",
                    DMDEDUC2 == 2 ~ "9-11th grade",
                    DMDEDUC2 == 3 ~ "High school graduate/GED",
                    DMDEDUC2 == 4 ~ "Some college or AA degree",
                    DMDEDUC2 == 5 ~ "College graduate or above",
                    TRUE ~ NA_character_
                ), levels = c("Less than 9th grade", "9-11th grade", "High school graduate/GED",
                             "Some college or AA degree", "College graduate or above")),
                Marital_Status = factor(dplyr::case_when(
                    DMDMARTZ == 1 ~ "Married",
                    DMDMARTZ == 2 ~ "Widowed",
                    DMDMARTZ == 3 ~ "Divorced",
                    DMDMARTZ == 4 ~ "Separated",
                    DMDMARTZ == 5 ~ "Never married",
                    DMDMARTZ == 6 ~ "Living with partner",
                    TRUE ~ NA_character_
                ), levels = c("Married", "Living with partner", "Never married", "Divorced", "Widowed", "Separated")),
                Race_Ethnicity = factor(dplyr::case_when(
                    RIDRETH1 == 1 ~ "Mexican American",
                    RIDRETH1 == 2 ~ "Other Hispanic",
                    RIDRETH1 == 3 ~ "Non-Hispanic White",
                    RIDRETH1 == 4 ~ "Non-Hispanic Black",
                    RIDRETH1 == 5 ~ "Other Race/Multi-Racial",
                    TRUE ~ NA_character_
                )),
                Smoking_Status = factor(dplyr::case_when(
                    SMQ020 == 2 ~ "Non-smoker",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
                    TRUE ~ NA_character_
                ), levels = c("Non-smoker", "Former smoker", "Current smoker")),
                Obesity = factor(dplyr::if_else(!is.na(BMXBMI) & BMXBMI >= 30, "Yes", "No", missing = "No"),
                                 levels = c("No", "Yes")),
                Diabetes_Mellitus = factor(dplyr::if_else(
                    DIQ010 == 1 | DIQ050 == 1 | DIQ070 == 1 |
                        (!is.na(LBXGH) & LBXGH >= 6.5) |
                        (!is.na(LBXSGL) & LBXSGL >= 126),
                    "Yes", "No", missing = "No"
                ), levels = c("No", "Yes")),
                Hypertension = factor(dplyr::if_else(BPQ020 == 1 | BPQ040A == 1, "Yes", "No", missing = "No"),
                                      levels = c("No", "Yes")),
                Chronic_Kidney_Disease = factor(dplyr::if_else((KIQ022 == 1 | KIQ025 == 1 | (!is.na(eGFR) & eGFR < 60)),
                                                               "Yes", "No", missing = "No"),
                                                levels = c("No", "Yes")),
                Cancer = factor(dplyr::if_else(MCQ220 == 1, "Yes", "No", missing = "No"),
                                levels = c("No", "Yes")),
                Cardiovascular_Disease = factor(dplyr::if_else(MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160D == 1,
                                                               "Yes", "No", missing = "No"),
                                                levels = c("No", "Yes"))
            )
    }
    
    rename_variables_safe <- function(data) {
        cat("Bắt đầu: Safe rename\n")
        rename_map <- c(
            ID = "SEQN",
            Age = "RIDAGEYR",
            Lymphocyte_percent = "LBXLYPCT",
            Neutrophil_percent = "LBXNEPCT",
            Albumin_g_L = "LBDSALSI"
        )
        safe_rename(data, rename_map)
    }
    
    create_all_features <- function(data) {
        cat("--- BẮT ĐẦU: FEATURE ENGINEERING ---\n")
        cols_to_check <- c("LBXSGL","SMQ040","DIQ050","DIQ070","BPQ040A","KIQ022","KIQ025","MCQ220",
                           "MCQ160B","MCQ160C","MCQ160D","MCQ160E","BPQ020",
                           "ALQ121","ALQ130","BMXBMI","BMXWAIST","BMXHIP","LBXGLU","LBDGLUSI",
                           "BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3","BPQ090D",
                           "LBDTRSI","LBDHDDSI","DMDEDUC2","DMDMARTZ","RIDRETH1","SMQ020")
        for (col in cols_to_check) if (!col %in% names(data)) data[[col]] <- NA
    
        data <- data %>%
            ensure_id() %>%
            derive_blood_indices() %>%
            define_diagnostic_criteria() %>%
            define_categorical_variables() %>%
            rename_variables_safe() %>%
            ensure_id()
    
        if (!("Age" %in% names(data))) data$Age <- as.numeric(data$RIDAGEYR)
        cat("--- HOÀN TẤT FEATURE ENGINEERING ---\n\n")
        data
    }
    
    # -------------------------
    # 3) OUTCOME TABLE FUNCTION (PSM)
    # -------------------------
    create_outcome_table <- function(data_pre, data_post, psm_vars, outcome_var,
                                     analysis_caption, control_level, treatment_level, Ns) {
    
        cat(paste("Tạo bảng outcome:", outcome_var, "\n"))
    
        unmatched_g0 <- data_pre %>% dplyr::filter(group == 0) %>% dplyr::pull(!!rlang::sym(outcome_var))
        unmatched_g1 <- data_pre %>% dplyr::filter(group == 1) %>% dplyr::pull(!!rlang::sym(outcome_var))
        val_g0_unmatched <- format_median_iqr(unmatched_g0)
        val_g1_unmatched <- format_median_iqr(unmatched_g1)
        p_unmatched <- tryCatch(stats::wilcox.test(unmatched_g0, unmatched_g1)$p.value, error = function(e) NA_real_)
    
        data_post_paired <- data_post %>% dplyr::arrange(subclass, group)
        matched_g0 <- data_post_paired %>% dplyr::filter(group == 0) %>% dplyr::pull(!!rlang::sym(outcome_var))
        matched_g1 <- data_post_paired %>% dplyr::filter(group == 1) %>% dplyr::pull(!!rlang::sym(outcome_var))
        val_g0_matched <- format_median_iqr(matched_g0)
        val_g1_matched <- format_median_iqr(matched_g1)
        p_matched <- tryCatch(stats::wilcox.test(matched_g0, matched_g1, paired = TRUE)$p.value, error = function(e) NA_real_)
    
        formula_adjusted <- stats::as.formula(paste(outcome_var, "~ group +", paste(psm_vars, collapse = " + ")))
        model_adj <- stats::lm(formula_adjusted, data = data_post, weights = data_post$weights)
        tidy_adj <- broom::tidy(model_adj) %>% dplyr::filter(term == "group")
        ci_adj <- tryCatch(stats::confint(model_adj), error = function(e) NULL)
    
        beta_formatted <- "—"
        p_adjusted <- NA_real_
        if (nrow(tidy_adj) == 1 && !is.null(ci_adj) && ("group" %in% rownames(ci_adj))) {
            beta_formatted <- paste0(
                format(round(tidy_adj$estimate, 2), nsmall = 2), " (",
                format(round(ci_adj["group", 1], 2), nsmall = 2), " - ",
                format(round(ci_adj["group", 2], 2), nsmall = 2), ")"
            )
            p_adjusted <- tidy_adj$p.value
        }
    
        outcome_df <- data.frame(
            Categories = c(paste0(outcome_var, ", Median (IQR)"), "Adjusted Beta (95% CI)"),
            Original_Control = c(val_g0_unmatched, "—"),
            Original_Treatment = c(val_g1_unmatched, "—"),
            Original_P = c(format_p_value(p_unmatched), "—"),
            Matched_Control = c(val_g0_matched, "—"),
            Matched_Treatment = c(val_g1_matched, beta_formatted),
            Matched_P = c(format_p_value(p_matched), format_p_value(p_adjusted)),
            stringsAsFactors = FALSE
        )
    
        flextable::flextable(outcome_df) %>%
            flextable::set_caption(analysis_caption) %>%
            flextable::set_header_labels(
                Categories = "Categories",
                Original_Control = paste0(control_level, " (N=", Ns$n_c_un, ")"),
                Original_Treatment = paste0(treatment_level, " (N=", Ns$n_t_un, ")"),
                Original_P = "P-value",
                Matched_Control = paste0(control_level, " (N=", Ns$n_c_m, ")"),
                Matched_Treatment = paste0(treatment_level, " (N=", Ns$n_t_m, ")"),
                Matched_P = "P-value"
            ) %>%
            flextable::add_header_row(values = c("", "Original Cohort", "Matched Cohort"),
                                      colwidths = c(1, 3, 3)) %>%
            flextable::theme_booktabs() %>%
            flextable::autofit() %>%
            flextable::align(j = 2:7, align = "center", part = "all") %>%
            flextable::align(j = 1, align = "left", part = "all") %>%
            fix_border_issues()
    }
    
    # -------------------------
    # 4) RUN PSM ANALYSIS
    # -------------------------
    run_psm_analysis <- function(data_input, group_var_string, treatment_level, control_level,
                                 psm_vars, outcome_var, output_filename, table_caption) {
    
        cat(paste0("\n--- BẮT ĐẦU PSM: ", group_var_string, " | ", treatment_level, " vs ", control_level, " ---\n"))
    
        check_required_cols(data_input, c("ID", group_var_string, psm_vars, outcome_var), context = "PSM inputs")
    
        data_pre <- data_input %>%
            dplyr::filter(!!rlang::sym(group_var_string) %in% c(treatment_level, control_level)) %>%
            dplyr::mutate(group = ifelse(!!rlang::sym(group_var_string) == treatment_level, 1, 0)) %>%
            dplyr::select(ID, group, dplyr::all_of(psm_vars), dplyr::all_of(outcome_var)) %>%
            tidyr::drop_na()
    
        Ns_list <- list(
            n_c_un = sum(data_pre$group == 0),
            n_t_un = sum(data_pre$group == 1)
        )
    
        if (Ns_list$n_c_un < 10 || Ns_list$n_t_un < 10) {
            cat("CẢNH BÁO: Không đủ quan sát. Dừng.\n")
            return(invisible(NULL))
        }
    
        formula_psm <- stats::as.formula(paste("group ~", paste(psm_vars, collapse = " + ")))
        psm_model <- MatchIt::matchit(formula = formula_psm, data = data_pre, method = "nearest",
                                      distance = "glm", caliper = 0.1, ratio = 1, replace = FALSE)
        data_post <- MatchIt::match.data(psm_model, data = data_pre)
    
        Ns_list$n_c_m <- sum(data_post$group == 0)
        Ns_list$n_t_m <- sum(data_post$group == 1)
    
        # --- Love Plot (Visual Check) ---
        plot_filename <- gsub("\\.docx$", ".png", output_filename)
        lp <- tryCatch({
            cobalt::love.plot(psm_model, stats = "m", binary = "std", abs = TRUE, thresholds = c(m = 0.2)) +
                ggplot2::theme_minimal(base_size = 12) +
                ggplot2::labs(title = paste("Covariate Balance:", treatment_level, "vs.", control_level),
                              x = "Absolute Standardized Mean Difference", y = NULL)
        }, error = function(e) NULL)
    
        if (!is.null(lp)) {
            ggplot2::ggsave(plot_filename, plot = lp, width = 8, height = max(6, length(psm_model$X) * 0.3), dpi = 300, bg = "white")
            cat("-> Đã lưu Love plot:", normalizePath(plot_filename), "\n")
        }
    
        # --- Word Output ---
        cat("Xuất Word...\n")
        doc <- officer::read_docx()
        ft_outcome <- create_outcome_table(data_pre, data_post, psm_vars, outcome_var, table_caption, control_level, treatment_level, Ns_list)
        doc <- flextable::body_add_flextable(doc, value = ft_outcome)
        base::print(doc, target = output_filename)
        cat("-> Đã lưu Word:", normalizePath(output_filename), "\n")
    
        invisible(list(psm_model = psm_model, data_pre = data_pre, data_post = data_post))
    }
    
    # -------------------------
    # 5) BMJ STYLE FOREST PLOT (3 MODELS)
    # -------------------------
    .logit_or_ci <- function(fit, term) {
        tt <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
        row <- tt[tt$term == term, , drop = FALSE]
        data.frame(OR = row$estimate, LCL = row$conf.low, UCL = row$conf.high, p = row$p.value)
    }
    
    .robust_or_ci <- function(glm_fit, term) {
        V <- sandwich::vcovHC(glm_fit, type = "HC0")
        ct <- lmtest::coeftest(glm_fit, vcov. = V)
        i <- which(rownames(ct) == term)
        beta <- ct[i, 1]; se <- ct[i, 2]; z <- stats::qnorm(0.975)
        data.frame(OR = exp(beta), LCL = exp(beta - z * se), UCL = exp(beta + z * se),
                   p = 2 * stats::pnorm(-abs(ct[i, 3])))
    }
    
    theme_bmj <- function(base_size = 11, base_family = "sans") {
        ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
            ggplot2::theme(
                axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
                axis.line.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(color = "black", face = "bold", size = 10, hjust = 0),
                axis.text.x = ggplot2::element_text(color = "black"),
                plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
                plot.subtitle = ggplot2::element_text(size = 10, hjust = 0),
                plot.caption = ggplot2::element_text(size = 8, hjust = 0),
                panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
                legend.position = "none"
            )
    }
    
    forest_three_models_bmj <- function(data, outcome, exposure, covars,
                                        labels = c("Multivariable Model", "PSM (Conditional Logit)", "IPTW"),
                                        xlim = NULL, show = TRUE, save_path = NULL, title = "") {
    
        cat("--- Calculating 3 Models for Forest Plot ---\n")
    
        # Model 1: Multivariable Logistic
        f_exp <- stats::as.formula(paste(outcome, "~", exposure, "+", paste(covars, collapse = " + ")))
        m1 <- stats::glm(f_exp, data = data, family = stats::binomial())
        r1 <- .logit_or_ci(m1, exposure)
        n1 <- stats::nobs(m1)
    
        # Model 2: PSM + Conditional Logistic
        f_ps <- stats::as.formula(paste(exposure, "~", paste(covars, collapse = " + ")))
        mi <- MatchIt::matchit(f_ps, data = data, method = "nearest", distance = "glm", ratio = 1, caliper = 0.1, replace = FALSE)
        dm <- MatchIt::match.data(mi, data = data)
        dm$subclass <- as.factor(dm$subclass)
        m2 <- survival::clogit(stats::as.formula(paste(outcome, "~", exposure, "+ strata(subclass)")), data = dm)
        tt2 <- broom::tidy(m2, conf.int = TRUE, exponentiate = TRUE)
        row2 <- tt2[tt2$term == exposure, , drop = FALSE]
        r2 <- data.frame(OR = row2$estimate, LCL = row2$conf.low, UCL = row2$conf.high, p = row2$p.value)
        n2 <- nrow(dm)
    
        # Model 3: IPTW
        ps_fit <- stats::glm(f_ps, data = data, family = stats::binomial())
        pscore <- stats::predict(ps_fit, type = "response")
        z <- stats::model.matrix(stats::as.formula(paste("~", exposure)), data = data)[, 2]
        p_t <- mean(z)
        w <- ifelse(z == 1, p_t / pscore, (1 - p_t) / (1 - pscore))
        f_out <- stats::as.formula(paste(outcome, "~", exposure))
        m3 <- stats::glm(f_out, data = data, family = stats::binomial(), weights = w)
        r3 <- .robust_or_ci(m3, exposure)
        n3 <- stats::nobs(m3)
    
        labels_with_n <- c(
            paste0(labels[1], "\n(N=", n1, ")"),
            paste0(labels[2], "\n(N=", n2, ")"),
            paste0(labels[3], "\n(N=", n3, ")")
        )
    
        res <- rbind(cbind(Model = labels_with_n[1], r1),
                     cbind(Model = labels_with_n[2], r2),
                     cbind(Model = labels_with_n[3], r3))
    
        res$Model <- factor(res$Model, levels = rev(labels_with_n))
        res$or_txt <- sprintf("%.2f (%.2f to %.2f)", res$OR, res$LCL, res$UCL)
        res$p_txt <- ifelse(res$p < 0.001, "<0.001", sprintf("%.3f", res$p))
    
        if (is.null(xlim)) {
            max_val <- max(res$UCL, na.rm = TRUE)
            min_val <- min(res$LCL, na.rm = TRUE)
            xlim <- c(max(0.1, min_val * 0.8), max_val * 1.5)
        }
    
        bmj_blue <- "#005689"
    
        range_x <- xlim[2] - xlim[1]
        pos_or <- xlim[2] + 0.05 * range_x
        pos_p <- xlim[2] + 0.5 * range_x
        final_xmax <- xlim[2] + 0.8 * range_x
    
        plot_breaks <- scales::breaks_extended(n = 5)(xlim)
    
        p <- ggplot2::ggplot(res, ggplot2::aes(y = Model, x = OR)) +
            ggplot2::geom_vline(xintercept = 1, linetype = "solid", color = "grey50", linewidth = 0.5) +
            ggplot2::geom_errorbarh(ggplot2::aes(xmin = LCL, xmax = UCL), height = 0.2, color = bmj_blue, linewidth = 0.8) +
            ggplot2::geom_point(size = 3.5, shape = 15, color = bmj_blue) +
            ggplot2::geom_text(ggplot2::aes(x = pos_or, label = or_txt), hjust = 0, size = 3.5, family = "sans") +
            ggplot2::geom_text(ggplot2::aes(x = pos_p, label = p_txt), hjust = 0, size = 3.5, family = "sans") +
            ggplot2::annotate("text", x = pos_or, y = length(labels) + 0.6, label = "OR (95% CI)", fontface = "bold", hjust = 0, size = 3.5) +
            ggplot2::annotate("text", x = pos_p, y = length(labels) + 0.6, label = "P-value", fontface = "bold", hjust = 0, size = 3.5) +
            ggplot2::scale_x_continuous(limits = c(xlim[1], final_xmax), breaks = plot_breaks) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::labs(title = title, x = "Odds Ratio (log scale)") +
            theme_bmj() +
            ggplot2::theme(plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10))
    
        if (!is.null(save_path)) {
            ggplot2::ggsave(save_path, p, width = 8, height = 3.5, dpi = 600)
            cat("-> Đã lưu Forest Plot (BMJ style):", normalizePath(save_path), "\n")
        }
        if (isTRUE(show)) print(p)
    }
    
    # ==============================================================================
    # MAIN EXECUTION
    # ==============================================================================
    
    # 1. Load Data
    if (!exists("data") || !is.data.frame(data)) {
        stop("LỖI: Biến 'data' chưa tồn tại. Vui lòng load dataset của bạn vào biến 'data' trước khi chạy.", call. = FALSE)
    }
    
    # 2. Pipeline
    data_cleaned <- apply_exclusion_criteria(data)
    data_featured <- create_all_features(data_cleaned)
    
    # 3. Define Variables
    psm_vars_shared <- c("Age","Gender","Race_Ethnicity","Marital_Status","Education_Level",
                         "Smoking_Status","Obesity","Diabetes_Mellitus","Hypertension",
                         "Chronic_Kidney_Disease","Cancer","Cardiovascular_Disease")
    
    cat("\n======================================================\n")
    cat("   PHẦN 1: PSM ANALYSES (TABLES + LOVE PLOTS)   \n")
    cat("======================================================\n")
    
    # --- A. MASLD vs Non-MASLD ---
    run_psm_analysis(
        data_input = data_featured,
        group_var_string = "masld_group",
        treatment_level = "MASLD",
        control_level = "non-MASLD",
        psm_vars = psm_vars_shared,
        outcome_var = "ANLR",
        output_filename = "Table_PSM_MASLD_vs_NonMASLD.docx",
        table_caption = "Comparison of ANLR between MASLD and non-MASLD"
    )
    
    # --- B. SLF vs Non-SLF (in MASLD) ---
    run_psm_analysis(
        data_input = data_featured %>% dplyr::filter(masld_group == "MASLD"),
        group_var_string = "slf_group",
        treatment_level = "SLF",
        control_level = "non-SLF",
        psm_vars = psm_vars_shared,
        outcome_var = "ANLR",
        output_filename = "Table_PSM_SLF_vs_NonSLF.docx",
        table_caption = "Comparison of ANLR between SLF and non-SLF (within MASLD)"
    )
    
    cat("\n======================================================\n")
    cat("   PHẦN 2: FOREST PLOTS (BMJ STYLE - 3 MODELS)   \n")
    cat("======================================================\n")
    
    # --- A. Forest Plot for MASLD (Exposure: High ANLR) ---
    df_forest_masld <- data_featured %>%
        dplyr::transmute(
            Outcome_Bin = as.integer(masld_group == "MASLD"),
            Exposure_Bin = as.integer(ANLR > stats::median(ANLR, na.rm = TRUE)),
            dplyr::across(dplyr::all_of(psm_vars_shared))
        ) %>%
        tidyr::drop_na()
    
    forest_three_models_bmj(
        data = df_forest_masld,
        outcome = "Outcome_Bin",
        exposure = "Exposure_Bin",
        covars = psm_vars_shared,
        xlim = c(0.7, 1.6),
        title = "Association between High ANLR and MASLD",
        save_path = "Forest_BMJ_MASLD.png"
    )
    
    # --- B. Forest Plot for SLF (Exposure: High ANLR) ---
    df_forest_slf <- data_featured %>%
        dplyr::filter(masld_group == "MASLD") %>%
        dplyr::transmute(
            Outcome_Bin = as.integer(slf_group == "SLF"),
            Exposure_Bin = as.integer(ANLR > stats::median(ANLR, na.rm = TRUE)),
            dplyr::across(dplyr::all_of(psm_vars_shared))
        ) %>%
        tidyr::drop_na()
    
    forest_three_models_bmj(
        data = df_forest_slf,
        outcome = "Outcome_Bin",
        exposure = "Exposure_Bin",
        covars = psm_vars_shared,
        xlim = c(0.5, 1),
        title = "Association between High ANLR and SLF (in MASLD)",
        save_path = "Forest_BMJ_SLF.png"
    )
    
    cat("\n--- ĐÃ HOÀN TẤT TOÀN BỘ QUY TRÌNH ---\n")
    
    ```