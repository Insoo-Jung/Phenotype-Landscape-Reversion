##### library #####

library(BoolNet)

##### Figure 3A #####

set.seed(123)

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