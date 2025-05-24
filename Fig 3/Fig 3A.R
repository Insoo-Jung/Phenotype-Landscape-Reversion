##### library #####

library(BoolNet)
library(ggplot2)

##### Figure 3A #####

# Change the file path #
network_folder = "C:/Users/Insoo Jung/Desktop/network_insoo/"

IO_Networks = list.files(network_folder, pattern = ".txt")
Error_Mat = c()
for(i in 1:length(IO_Networks)){
  print(i)
  net = loadNetwork(paste0(network_folder,IO_Networks[i]))
  L_F = Large_L_F(net)
  L_B = Large_L_B(net)
  number = min(10000, 2^length(net$genes))
  set.seed(123)
  at = Large_Attr(net, number)
  for(j in 1:20){
    error = 0
    for(k in 1:dim(at$Attr_Mat)[2]){
      vector = at$Attr_Mat[,k]
      result = MM(L_F, j) %*% MM(L_B, j) %*% matrix(vector, ncol = 1)
      error = error + (sum(abs(result - vector))/(2 * length(net$genes))) * at$Attr_Size[k]
    }
    Error_Mat = c(Error_Mat, error / (sum(at$Attr_Size)))
  }
}

Network_Info = c()
for(i in 1:length(IO_Networks)){
  Network_Info = c(Network_Info, rep(paste0(strsplit(IO_Networks[[i]], split = ".txt")[[1]][1]), 20))
}
Steps = as.character(rep(seq(1,20,1), length(IO_Networks)))
large_Error = data.frame(Error_Mat, Steps, Network_Info)

save(large_Error, file = "Fig_3A_data.RData")

Fig_3A = ggplot(large_Error, aes(x = Steps, y = 1 - Error_Mat)) + geom_boxplot(notch = TRUE, fill = "#3399FF", outlier.shape = NA) + ylab("Accuracy") + theme_bw() + theme(legend.position = 'none', axis.text.x = element_text(size = 30, hjust = 1), axis.text.y = element_text(size = 30), plot.title = element_text(size = 30, hjust = 0.5), axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + xlab("Transition Steps") + scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")) + geom_jitter(width = 0.1, alpha = 0.3)