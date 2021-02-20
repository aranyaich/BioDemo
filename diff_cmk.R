

# 01 load data
raw_path = "C:/Users/alan/Documents/WeChat Files/wxid_b53hwbekwcxb22/FileStorage/File/2020-12/Õª±‰“¿¿µ–‘-CMK.csv"
raw_data = read.csv(raw_path)


# 02 clean data
exp_data <- raw_data[,7:dim(raw_data)[2]]

group_data = factor(raw_data[,"TERT.Mutation"])
group1_index = grep("WT",group_data)
group2_index = grep("C228T",group_data)


# 03 count the logfc and p-value

# logfc<-function(x){
#   log(mean(x[group1_index])/mean(x[group2_index]))
# }
logfc <- function(x){
  mean(x[group1_index])-mean(x[group2_index])
}
wil_p <- function(x){
  wilcox.test(x[group1_index],x[group2_index])$p.value
}

test_data = t(exp_data)
gene_fc = apply(test_data,1,logfc) # count log fold change for each gene
gene_pvalue = apply(test_data,1,wil_p) # count the p value for each gene
  

# 04 draw plot

library("ggplot2")

get_class <- function(x){
  if (!is.na(x["logfc"])&&x["pval"]<0.05) {
    if(x["logfc"]>0){
      return("up")
      }else{return("down")}
  }else{
    return("flat")
  }
}
combined_data = data.frame(logfc = gene_fc,pval = gene_pvalue)
combined_data$logfc[is.na(combined_data$logfc)]<-0
write.csv(combined_data,"C:/Users/alan/Documents/WeChat Files/wxid_b53hwbekwcxb22/FileStorage/File/2020-12/results_cmk.csv")

combined_data$order_fc = -rank(combined_data$logfc)
combined_data$class = apply(combined_data,1,get_class)


plot_data<- combined_data
ggplot(plot_data, aes(x=order_fc, y=logfc)) +  
  geom_point(shape=1,aes(color=class)) +
  labs(x="genes",y="log2 fold change") + theme_bw()


