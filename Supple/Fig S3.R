##### library #####

library(BoolNet)

##### Figure S3 #####

# Change the file path #
network_folder = "C:/Users/Insoo Jung/Desktop/network_insoo/"

set.seed(123)

IO_Networks = list.files(network_folder, pattern = ".txt")
O_net = c()
for(i in 1:length(IO_Networks)){
  net = loadNetwork(paste0(network_folder,IO_Networks[i]))
  Adj = matrix(0, length(net$genes), length(net$genes))
  for(j in 1:length(net$genes)){
    inputs = net$interactions[[j]]$input
    Adj[j,inputs] = 1
  }
  if(length(which(colSums(Adj) == 0)) > 0){
    O_net = c(O_net, IO_Networks[i])
  }
}

Time_Cost_Algo = c()
Time_Cost_Simul = c()
for(i in 1:length(O_net)){
  net = loadNetwork(paste0(network_folder, O_net[i]))
  N_Num = length(net$genes)
  print(i)
  
  at_start_time = Sys.time()
  attr = Large_Attr(net,min(10000,(2^N_Num)))
  at_end_time = Sys.time()
  
  Large_ARC_Search(net, attr$Attr_Mat, attr$Attr_Size, 1, 0, seq(1,N_Num,1))
  algo_end_time = Sys.time()
  
  Time_Cost_Algo = c(Time_Cost_Algo, difftime(algo_end_time, at_start_time, units = "mins"))
  Time_Cost_Simul = c(Time_Cost_Simul, 2 * N_Num * difftime(at_end_time, at_start_time, units = "mins"))
}

Ratio = Time_Cost_Simul / Time_Cost_Algo
ratio = Ratio[c(22,25,4,32,20,10,17,2,3,33,30,11,19,7,23,12,29,6,9,31,15,21,1,16,27,13,24,28,18,5,14,8,26)]
Index = seq(1,33,1)
S4_Data = data.frame(ratio, Index)
save(S4_Data, file = "Fig_S3_data.RData")

ggplot(S4_Data, aes(x = Index, y = ratio)) + geom_bar(stat = 'identity') + theme_bw() + ylab("Expected Time Ratio") + xlab("Network Index")
