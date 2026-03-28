library(readxl)
library(accucor)
library(openxlsx)
library(tidyverse)
rm(list = ls())

# 定义函数,质荷比-强度匹配，可更改ppm阈值,默认为5ppm
Extractint <- function(Intensity, Intensitys, next_Intensitys, thr=5*10^(-6)){
  results <- abs((Intensitys-Intensity)/Intensity)
  potential_Intensitys <- next_Intensitys[results < thr]
  if(length(potential_Intensitys) == 0){
    return(NA)
  }
  return(max(potential_Intensitys))
}
# 定义函数,过滤成像效果不好的离子，可更改阈值
Scanint <- function(df, Intensity=10, num=2){
  remove_ids <- NULL
  
  select_df <- df[df$IsotopeLabel == "C12 PARENT", ]
  for(row in 1:nrow(select_df)){
    if(sum(select_df[row, c("Unlabel_01", "Unlabel_02", "Unlabel_03", "Glc_Label_01", "Glc_Label_02", "Glc_Label_03")] < Intensity) > num){
      remove_ids <- c(remove_ids, unlist(select_df[row, "ID"]))
    }
  }
  
  remain_df <- df[df$ID %in% remove_ids == F, ]
  return(remain_df)
}
# 定义函数,保留MSI可见的标记离子
Filterint <- function(df, Intensity=10){
  new_df <- df[, 1:6]
  for(col in c("Unlabel_01", "Unlabel_02", "Unlabel_03", "Glc_Label_01", "Glc_Label_02", "Glc_Label_03")){
    Intensitys <- unlist(df[, col])
    new_df[, col] <- ifelse(Intensitys == 0, 1, Intensitys)
  }
  
  remove_ids2 <- NULL
  remove_ids3 <- NULL
  
  select_df <- new_df[new_df$IsotopeLabel == "C13-label-2",]
  for(row in 1:nrow(select_df)){
    sum_Intensity1 <- sum(unlist(select_df[row, c("Glc_Label_01", "Glc_Label_02", "Glc_Label_03")]))
    sum_Intensity2 <- sum(unlist(select_df[row, c("Unlabel_01", "Unlabel_02", "Unlabel_03")]))
    if(unlist(select_df[row, "ID"]) == 2){
      print(sum_Intensity1)
      print(sum_Intensity2)
      print(sum_Intensity1/sum_Intensity2)
    }
    if(sum_Intensity1/sum_Intensity2 < Intensity){
      remove_ids2 <- c(remove_ids2, unlist(select_df[row, "ID"]))
    }
  }
  
  select_df <- new_df[new_df$IsotopeLabel == "C13-label-3",]
  for(row in 1:nrow(select_df)){
    sum_Intensity1 <- sum(unlist(select_df[row, c("Glc_Label_01", "Glc_Label_02", "Glc_Label_03")]))
    sum_Intensity2 <- sum(unlist(select_df[row, c("Unlabel_01", "Unlabel_02", "Unlabel_03")]))
    if(sum_Intensity1/sum_Intensity2 < Intensity){
      remove_ids3 <- c(remove_ids3, unlist(select_df[row, "ID"]))
    }
  }
  
  remove_ids <- NULL
  for(id in remove_ids2){
    if(id %in% remove_ids3){
      remove_ids <- c(remove_ids, id)
    }
  }
  print(remove_ids)
  remain_df <- df[(df$ID %in% remove_ids) == F,]
  return(remain_df)
}

# 导入MSI数据及Target列表
# MSI数据更改名称为"Unlabel_01", "Unlabel_02", "Unlabel_03", "Glc_Label_01", "Glc_Label_02", "Glc_Label_03"
MSIdata <- read_xlsx("mydata.xlsx", sheet ="intensity")
Target <- read_xlsx("mydata.xlsx", sheet ="target")

# 强度提取匹配
for(i in seq(1, ncol(MSIdata), 2)){
  print(i)
  col <- colnames(MSIdata)[i]
  select_Intensitys <- apply(Target[, 1, drop=F], 1, Extractint, unlist(MSIdata[, i]), unlist(MSIdata[, i+1]))
  Target[, col] <- select_Intensitys
  Target[is.na(Target[, col]), col] <- 0
}
Target1 <- Scanint(Target,Intensity=100, num=2)
Target2 <- Filterint(Target1, Intensity=10)
# 离子强度匹配及筛选结果导出，为“Result_Int”
Compound_new <- Target2 %>% select(Compound,ID) %>% mutate(Compound_new = pmap_chr(., str_c, sep = "_")) %>% select(Compound_new)
Target3 <- cbind(Compound_new, Target2 %>% select(Formula,IsotopeLabel,Unlabel_01:Glc_Label_03)) %>% rename(Compound = Compound_new)
Result_Int = list("Extract_Int" = Target,"Scan_Int" = Target1,"Filter_Int"=Target2,"Input_Int"=Target3)
write.xlsx(Result_Int,file = "Result.xlsx")
# 天然同位素校正，结果自动导出为“Result_Int_corrected”
Corrected <- natural_abundance_correction(path = "Result.xlsx",sheet = "Input_Int",resolution = 140000,purity = 0.99)



