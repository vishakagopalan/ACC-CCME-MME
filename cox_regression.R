correct_feature_names <- function(feature_df,feature_names) {
    new_feature_name_vec <- colnames(feature_df)
    names(new_feature_name_vec) <- colnames(feature_df)
    for (feature in colnames(feature_df)) {
        if (feature %in% feature_names) {
            new_name <- gsub("-|\\ |:|/","_",feature)
            new_feature_name_vec[feature] <- new_name
        } 
    }

    colnames(feature_df) <- new_feature_name_vec
    new_feature_name_vec <- new_feature_name_vec[feature_names]
    return(list("feature_df"=feature_df,"new_names"=new_feature_name_vec))
}

cox_regression <- function( sample_survival_df, sample_features_df, survival_data_id_col="sample_id",
feature_data_id_col="sample_id", coef_suffix="",
model_covariates=c("strata(gender)","age_at_index"),
status_str="status", duration_str="days_to_event",
low_high_percentiles=c(0.25,0.75),km_features_to_plot=NULL,do_forest_plot=F ) {
    features <- setdiff(colnames(sample_features_df),feature_data_id_col)

    out <- correct_feature_names(sample_features_df,features)

    if (survival_data_id_col == feature_data_id_col) {
        cox_df <- merge( sample_survival_df, out$feature_df, by=feature_data_id_col )
    } else {
        cox_df <- merge( sample_survival_df, out$feature_df, by.x=survival_data_id_col, by.y=feature_data_id_col )
    }

    regression_str_lhs <- paste0("Surv(",duration_str,",",status_str,")")
    p_value_df <- data.frame()

    if (!is.null(model_covariates)) {
        for (idx in 1:length(model_covariates)) {
            model_covariates[idx] <- paste0("`",model_covariates[idx],"`")
        }
    }
    idx <- 1
    forest_plot_list <- list()
    for (feature in out$new_names) {
        print(paste(idx,length(out$new_names),sep="/"))
        flush.console()
        idx <- idx + 1
        if (sum(is.na(cox_df[[feature]])) == nrow(cox_df))
            next
        
        regression_str_rhs <- feature
        if (!is.null(model_covariates)) {
            regression_str_rhs <- paste(c(regression_str_rhs,model_covariates),collapse="+")
        }

        regression_str <- paste(regression_str_lhs,regression_str_rhs,sep="~")
        other_features <- setdiff(out$new_names,feature)
        sub_df <- dplyr::select(cox_df,-all_of(other_features))

        print(regression_str)
        cox_obj <- coxph(as.formula(regression_str), data=sub_df)
        if (do_forest_plot)
            forest_plot_list[[feature]] <- list("cox_obj"=cox_obj,"for_cox_df"=sub_df )
        cox_summ <- summary(cox_obj)
        p_value <- cox_summ$coefficients[paste0(feature,coef_suffix),"Pr(>|z|)"]
        hazard_ratio <- cox_summ$coefficients[paste0(feature,coef_suffix),"exp(coef)"]
        logrank_p_value <- cox_summ$logtest[3]
        
        p_value_df <- rbind(p_value_df,data.frame(feature=feature,
                                                  hr_p_value=p_value,hazard_ratio=hazard_ratio,
                                                  logrank_p_value=logrank_p_value))

    }    

    p_value_df <- merge(p_value_df, tibble::enframe(out$new_names,name="actual_feature",value="renamed_feature"),
    by.x="feature",by.y="renamed_feature") %>% dplyr::select(-feature) %>% dplyr::rename(feature=actual_feature)

    to_return_list <- list("regression_df"=p_value_df)
    if (do_forest_plot)
        to_return_list[["forest_info"]] <- forest_plot_list

    if (!is.null(km_features_to_plot)) {
        plot_list <- list()
        covariate_cols <- gsub("strata(\\(..*?\\))","\\1",model_covariates) %>% gsub("\\(|\\)","",.)
        km_df <- dplyr::select(cox_df,all_of(c(duration_str,status_str,covariate_cols))) %>% as.data.frame
        for (km_feature in km_features_to_plot) {
            column <- paste(km_feature,"status",sep="_")
            print(column)
            flush.console()
            regression_str_rhs <- paste(c(km_feature,model_covariates),collapse="+")
            regression_str <- paste(regression_str_lhs,regression_str_rhs,sep="~")
            cox_df <- mutate(cox_df,!!sym(column):=case_when(!!sym(km_feature) > quantile(!!sym(km_feature),low_high_percentiles[2]) ~ "High",
                                                                             !!sym(km_feature) < quantile(!!sym(km_feature),low_high_percentiles[1]) ~ "Low",
                                                                             TRUE ~ "Medium"))
            km_df[[km_feature]] <- cox_df[[column]]
        }

        to_return_list[["km_df"]]=km_df
    }

    return(to_return_list)

}
