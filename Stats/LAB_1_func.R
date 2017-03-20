difference = matrix(0, 400,1)
i=1
for (n in 100:1000){ 
  y = rnorm(n)
  y_bar = mean(y)

  
  s = sum((y - y_bar)**2)/(n-1)
  s_bar = sum((y - y_bar)**2)/(n)
  
  difference[i] = abs(s - s_bar)
  i=i+1

}



plot(1:length(s_keep), difference)
