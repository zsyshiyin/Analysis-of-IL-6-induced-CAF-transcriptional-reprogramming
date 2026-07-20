# =============================================================================
# 生信分析脚本：ssGSEA评分与临床数据关联分析（含评分缩放版）
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
work_dir <- "D:/zsy/SX/Fomal-final/13-TCGA-clinical"
setwd(work_dir)

output_dir <- "survival_analysis-TNM"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

clinical_file <- "TCGA_clinical_processed.csv"
score_file <- "ssGSEA_scores_diagnostic.csv"

clinical_id_col <- "Id"
clinical_futime <- "futime"
clinical_fustat <- "fustat"
clinical_vars <- c("age", "gender", "grade", "stage", "T", "M", "N")

score_id_col <- "Sample"
score_value_col <- "ssGSEA_Score"

# 评分预处理设置
score_transform <- "scale_10"  # "none", "zscore", "scale_10"

score_split_method <- "median"

km_width <- 4.5
km_height <- 4.5
km_title <- "Survival Analysis by ssGSEA Score"
km_xlab <- "Time (days)"
km_ylab <- "Survival Probability"
km_palette <- c("#25377F", "#C60036")
km_font_base <- 18
km_font_title <- 22
km_risk_table <- TRUE

cox_vars <- c("age", "gender", "grade", "stage", "T", "M", "N")
cox_show_univariate <- TRUE
cox_show_multivariate <- TRUE

output_km <- file.path(output_dir, "KM_curve.pdf")
output_km_with_cox <- file.path(output_dir, "KM_curve_with_cox.pdf")
output_cox <- file.path(output_dir, "COX_regression_results.csv")
output_data <- file.path(output_dir, "merged_survival_data.csv")
output_summary <- file.path(output_dir, "survival_analysis_summary.txt")

# ====================== 加载必要的包 ======================
cat("\n========== 加载必要的R包 ==========\n")

packages <- c("survival", "survminer", "ggplot2", "dplyr", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取数据 ======================
cat("\n========== 读取数据 ==========\n")

clinical_data <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
score_data <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)

# ====================== ID标准化和匹配 ======================
cat("\n========== ID标准化和匹配 ==========\n")

clinical_data$PatientID <- substr(clinical_data[[clinical_id_col]], 1, 12)
score_data$PatientID <- substr(score_data[[score_id_col]], 1, 12)

common_patients <- intersect(clinical_data$PatientID, score_data$PatientID)
cat(sprintf("\n共同患者数量: %d\n", length(common_patients)))

if (length(common_patients) == 0) {
  stop("错误：没有匹配的患者ID！")
}

merged_data <- clinical_data %>%
  filter(PatientID %in% common_patients) %>%
  left_join(score_data[, c("PatientID", score_value_col)], by = "PatientID") %>%
  distinct(PatientID, .keep_all = TRUE)

# ====================== 数据预处理和评分缩放 ======================
cat("\n========== 数据预处理和评分缩放 ==========\n")

colnames(merged_data)[colnames(merged_data) == score_value_col] <- "ssGSEA_Score_Original"

merged_data$futime <- as.numeric(merged_data[[clinical_futime]])
merged_data$fustat <- as.numeric(merged_data[[clinical_fustat]])

valid_survival <- merged_data %>%
  filter(!is.na(futime) & !is.na(fustat) & futime > 0)

cat(sprintf("\n有效生存数据样本数: %d\n", nrow(valid_survival)))

if (nrow(valid_survival) == 0) {
  stop("错误：没有有效的生存数据！")
}

# 显示原始评分统计
cat("\n原始ssGSEA评分统计:\n")
cat(sprintf("  范围: [%.6f, %.6f]\n", 
            min(valid_survival$ssGSEA_Score_Original), 
            max(valid_survival$ssGSEA_Score_Original)))
cat(sprintf("  均值: %.6f\n", mean(valid_survival$ssGSEA_Score_Original)))
cat(sprintf("  标准差: %.6f\n", sd(valid_survival$ssGSEA_Score_Original)))

# 评分转换
if (score_transform == "zscore") {
  valid_survival$ssGSEA_Score <- scale(valid_survival$ssGSEA_Score_Original)
  valid_survival$ssGSEA_Score <- as.numeric(valid_survival$ssGSEA_Score)
  cat("\n已应用Z-score标准化\n")
  cat("转换后评分统计:\n")
  cat(sprintf("  范围: [%.4f, %.4f]\n", 
              min(valid_survival$ssGSEA_Score), 
              max(valid_survival$ssGSEA_Score)))
  cat(sprintf("  均值: %.4f\n", mean(valid_survival$ssGSEA_Score)))
  cat(sprintf("  标准差: %.4f\n", sd(valid_survival$ssGSEA_Score)))
  
} else if (score_transform == "scale_10") {
  valid_survival$ssGSEA_Score <- valid_survival$ssGSEA_Score_Original * 10
  cat("\n已应用乘以10缩放\n")
  cat("转换后评分统计:\n")
  cat(sprintf("  范围: [%.4f, %.4f]\n", 
              min(valid_survival$ssGSEA_Score), 
              max(valid_survival$ssGSEA_Score)))
  cat(sprintf("  均值: %.4f\n", mean(valid_survival$ssGSEA_Score)))
  cat(sprintf("  标准差: %.4f\n", sd(valid_survival$ssGSEA_Score)))
  
} else {
  valid_survival$ssGSEA_Score <- valid_survival$ssGSEA_Score_Original
  cat("\n使用原始评分\n")
}

# 移除原始评分列
valid_survival$ssGSEA_Score_Original <- NULL

# ====================== 评分分组 ======================
cat("\n========== 评分分组 ==================\n")

if (score_split_method == "median") {
  score_cutoff <- median(valid_survival$ssGSEA_Score, na.rm = TRUE)
  valid_survival$Score_group <- ifelse(valid_survival$ssGSEA_Score >= score_cutoff, "High", "Low")
  valid_survival$Score_group <- factor(valid_survival$Score_group, levels = c("Low", "High"))
  cat(sprintf("使用中位数分组 (cutoff = %.4f)\n", score_cutoff))
}

cat("\n评分分组统计:\n")
print(table(valid_survival$Score_group))

# ====================== KM生存曲线 ======================
cat("\n========== 绘制KM生存曲线 ==========\n")

surv_obj <- Surv(time = valid_survival$futime, event = valid_survival$fustat)
km_fit <- survfit(surv_obj ~ Score_group, data = valid_survival)

logrank_test <- survdiff(surv_obj ~ Score_group, data = valid_survival)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = 1)
cat(sprintf("\nLog-rank检验 p值: %.4e\n", logrank_p))

# 计算中位生存时间
km_summary <- summary(km_fit)
median_surv <- km_summary$table[, "median"]
cat("\n中位生存时间:\n")
print(median_surv)

# 绘制KM曲线
km_plot <- ggsurvplot(km_fit,
                      data = valid_survival,
                      pval = TRUE,
                      pval.method = TRUE,
                      conf.int = TRUE,
                      risk.table = km_risk_table,
                      risk.table.col = "strata",
                      palette = km_palette,
                      xlab = km_xlab,
                      ylab = km_ylab,
                      title = km_title,
                      legend.title = "Score Group",
                      legend.labs = levels(valid_survival$Score_group),
                      ggtheme = theme_minimal())

pdf(output_km, width = km_width, height = km_height)
print(km_plot)
dev.off()
cat("KM曲线已保存到:", output_km, "\n")

# ====================== COX回归分析 ======================
cat("\n========== COX回归分析 ==========\n")

cox_data <- valid_survival %>%
  select(futime, fustat, ssGSEA_Score, all_of(cox_vars)) %>%
  na.omit()

cat(sprintf("COX回归可用样本数: %d\n", nrow(cox_data)))

# 处理分类变量
prepare_categorical <- function(data, var) {
  if (var %in% colnames(data)) {
    data[[var]] <- as.factor(data[[var]])
    levels_count <- length(unique(na.omit(data[[var]])))
    if (levels_count < 2) {
      return(NULL)
    }
    return(data[[var]])
  }
  return(NULL)
}

# 单因素COX回归
if (cox_show_univariate) {
  cat("\n单因素COX回归分析:\n")
  univ_results <- data.frame()
  
  for (var in c("ssGSEA_Score", cox_vars)) {
    if (var %in% colnames(cox_data)) {
      tryCatch({
        formula <- as.formula(paste("Surv(futime, fustat) ~", var))
        fit <- coxph(formula, data = cox_data)
        summary_fit <- summary(fit)
        coef <- summary_fit$coefficients
        
        for (i in 1:nrow(coef)) {
          var_name <- rownames(coef)[i]
          conf_int <- tryCatch(exp(confint(fit))[i, ], error = function(e) c(NA, NA))
          
          univ_results <- rbind(univ_results, data.frame(
            Variable = var_name,
            HR = coef[i, "exp(coef)"],
            HR_lower = conf_int[1],
            HR_upper = conf_int[2],
            P_value = coef[i, "Pr(>|z|)"],
            stringsAsFactors = FALSE
          ))
          
          hr_note <- ""
          if (var == "ssGSEA_Score") {
            if (score_transform == "scale_10") {
              hr_note <- " (原始评分每增加0.1的风险比)"
            } else if (score_transform == "zscore") {
              hr_note <- " (每增加1个标准差的风险比)"
            }
          }
          
          cat(sprintf("  %s%s: HR = %.3f (%.3f-%.3f), p = %.4e\n", 
                      var_name, hr_note, coef[i, "exp(coef)"], 
                      conf_int[1], conf_int[2],
                      coef[i, "Pr(>|z|)"]))
        }
      }, error = function(e) {
        cat(sprintf("  变量 %s 分析失败: %s\n", var, e$message))
      })
    }
  }
}

# 多因素COX回归
if (cox_show_multivariate) {
  cat("\n多因素COX回归分析:\n")
  
  vars_to_include <- c("ssGSEA_Score")
  for (var in cox_vars) {
    if (var %in% colnames(cox_data)) {
      if (var %in% c("gender", "grade", "stage", "T", "M", "N")) {
        temp_var <- prepare_categorical(cox_data, var)
        if (!is.null(temp_var)) {
          vars_to_include <- c(vars_to_include, var)
        }
      } else {
        vars_to_include <- c(vars_to_include, var)
      }
    }
  }
  
  if (length(vars_to_include) > 1) {
    formula_multi <- as.formula(paste("Surv(futime, fustat) ~", 
                                      paste(vars_to_include, collapse = " + ")))
    
    tryCatch({
      multi_fit <- coxph(formula_multi, data = cox_data)
      summary_multi <- summary(multi_fit)
      
      multi_results <- data.frame()
      coef_multi <- summary_multi$coefficients
      
      for (i in 1:nrow(coef_multi)) {
        conf_int_multi <- tryCatch(exp(confint(multi_fit))[i, ], error = function(e) c(NA, NA))
        
        multi_results <- rbind(multi_results, data.frame(
          Variable = rownames(coef_multi)[i],
          HR = coef_multi[i, "exp(coef)"],
          HR_lower = conf_int_multi[1],
          HR_upper = conf_int_multi[2],
          P_value = coef_multi[i, "Pr(>|z|)"],
          stringsAsFactors = FALSE
        ))
        
        hr_note <- ""
        if (grepl("ssGSEA_Score", rownames(coef_multi)[i])) {
          if (score_transform == "scale_10") {
            hr_note <- " (原始评分每增加0.1的风险比)"
          } else if (score_transform == "zscore") {
            hr_note <- " (每增加1个标准差的风险比)"
          }
        }
        
        cat(sprintf("  %s%s: HR = %.3f (%.3f-%.3f), p = %.4e\n", 
                    rownames(coef_multi)[i], hr_note,
                    coef_multi[i, "exp(coef)"],
                    conf_int_multi[1],
                    conf_int_multi[2],
                    coef_multi[i, "Pr(>|z|)"]))
      }
    }, error = function(e) {
      cat("多因素COX回归失败:", e$message, "\n")
      multi_results <- data.frame()
    })
  } else {
    cat("没有足够的变量进行多因素COX回归\n")
    multi_results <- data.frame()
  }
}

# ====================== 保存COX回归结果 ======================
cat("\n========== 保存COX回归结果 ==========\n")

all_cox_results <- data.frame()

if (exists("univ_results") && nrow(univ_results) > 0) {
  all_cox_results <- rbind(all_cox_results,
                           data.frame(Type = "Univariate", univ_results))
}

if (exists("multi_results") && nrow(multi_results) > 0) {
  all_cox_results <- rbind(all_cox_results,
                           data.frame(Type = "Multivariate", multi_results))
}

if (nrow(all_cox_results) > 0) {
  write.csv(all_cox_results, file = output_cox, row.names = FALSE)
  cat("COX回归结果已保存到:", output_cox, "\n")
} else {
  cat("警告：没有生成COX回归结果！\n")
}

# ====================== 在KM曲线上添加COX结果 ======================
cat("\n========== 生成带COX结果的KM曲线 ==========\n")

cox_score_result <- subset(all_cox_results, 
                           grepl("ssGSEA_Score", Variable) & Type == "Multivariate")
if (nrow(cox_score_result) == 0) {
  cox_score_result <- subset(all_cox_results, 
                             grepl("ssGSEA_Score", Variable) & Type == "Univariate")
}

if (nrow(cox_score_result) > 0) {
  hr_note <- ""
  if (score_transform == "scale_10") {
    hr_note <- " (per 0.1 increase)"
  } else if (score_transform == "zscore") {
    hr_note <- " (per SD increase)"
  }
  
  cox_text <- sprintf("COX%s: HR = %.2f (%.2f-%.2f), p = %.3e", 
                      hr_note,
                      cox_score_result$HR[1],
                      cox_score_result$HR_lower[1],
                      cox_score_result$HR_upper[1],
                      cox_score_result$P_value[1])
} else {
  cox_text <- sprintf("Log-rank p = %.3e", logrank_p)
}

km_plot_with_cox <- km_plot
km_plot_with_cox$plot <- km_plot_with_cox$plot + 
  annotate("text", x = Inf, y = 0.2, 
           label = cox_text,
           hjust = 1.1, size = 3, color = "darkred")

output_km_with_cox <- file.path(output_dir, "KM_curve_with_cox.pdf")
pdf(output_km_with_cox, width = km_width, height = km_height)
print(km_plot_with_cox)
dev.off()

cat("带COX结果的KM曲线已保存到:", output_km_with_cox, "\n")

# ====================== 保存合并后的数据 ======================
cat("\n========== 保存合并数据 ==========\n")

write.csv(valid_survival, file = output_data, row.names = FALSE)
cat("合并后的数据已保存到:", output_data, "\n")

# ====================== 生成分析总结 ======================
cat("\n========== 生成分析总结 ==========\n")

sink(output_summary)

cat("========================================\n")
cat("      ssGSEA评分生存分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、数据概况\n")
cat("------------\n")
cat(sprintf("总样本数: %d\n", nrow(valid_survival)))
cat(sprintf("事件数: %d\n", sum(valid_survival$fustat == 1)))
cat(sprintf("删失数: %d\n", sum(valid_survival$fustat == 0)))
cat(sprintf("中位随访时间: %.1f 天\n", median(valid_survival$futime)))

cat("\n二、评分预处理\n")
cat("--------------\n")
cat(sprintf("转换方法: %s\n", score_transform))
if (score_transform == "scale_10") {
  cat("说明: 评分乘以10，使HR表示原始评分每增加0.1的风险比\n")
} else if (score_transform == "zscore") {
  cat("说明: Z-score标准化，使HR表示每增加1个标准差的风险比\n")
} else {
  cat("说明: 使用原始评分（范围-1到1）\n")
}

cat("\n三、评分分组\n")
cat("------------\n")
cat(sprintf("分组方法: %s\n", score_split_method))
if (exists("score_cutoff")) {
  cat(sprintf("阈值: %.4f\n", score_cutoff))
}
print(table(valid_survival$Score_group))

cat("\n四、生存分析结果\n")
cat("----------------\n")
cat(sprintf("Log-rank检验 p值: %.4e\n", logrank_p))
cat("\n中位生存时间:\n")
print(median_surv)

if (exists("cox_score_result") && nrow(cox_score_result) > 0) {
  cat("\n五、COX回归结果（ssGSEA评分）\n")
  cat("-------------------------------\n")
  if (score_transform == "scale_10") {
    cat(sprintf("HR (每增加0.1): %.3f (%.3f-%.3f)\n", 
                cox_score_result$HR[1],
                cox_score_result$HR_lower[1],
                cox_score_result$HR_upper[1]))
  } else if (score_transform == "zscore") {
    cat(sprintf("HR (每增加1个标准差): %.3f (%.3f-%.3f)\n", 
                cox_score_result$HR[1],
                cox_score_result$HR_lower[1],
                cox_score_result$HR_upper[1]))
  } else {
    cat(sprintf("HR: %.3f (%.3f-%.3f)\n", 
                cox_score_result$HR[1],
                cox_score_result$HR_lower[1],
                cox_score_result$HR_upper[1]))
  }
  cat(sprintf("p值: %.4e\n", cox_score_result$P_value[1]))
}

cat("\n六、输出文件列表\n")
cat("----------------\n")
cat("1.", basename(output_km), "- KM生存曲线\n")
cat("2.", basename(output_km_with_cox), "- 带COX结果的KM曲线\n")
cat("3.", basename(output_cox), "- COX回归结果\n")
cat("4.", basename(output_data), "- 合并后的数据\n")

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("分析总结已保存到:", output_summary, "\n")
cat("\n=================== 分析完成 ===================\n")