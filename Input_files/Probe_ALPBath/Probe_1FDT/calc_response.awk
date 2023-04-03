
BEGIN{
  COR=1000
  DT=0.01
  
  tot_files = 0
}
{
  sign = 1
  if (FILENAME=="colloid.dat") sign = -1
  
  response[(FNR-1)%(COR+1)]+= sign*$4
  
  if ((FNR-1)%(COR+1) == 0 ) tot_files ++
}
END{
  for (t=0; t<COR; t++) {
    if (t>0) print t*DT, response[t]/((double)tot_files)*2.0
    else print t*DT, -0.5
    
  }
  
  
}

