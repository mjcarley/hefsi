function [x,t]=readtri(file)

  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%d", 2) ;

  nn = dat(1) ; ne = dat(2) ;
  
  x = fscanf(fid, "%f", nn*5) ;
  x = reshape(x, 5, nn)' ;

  t = fscanf(fid, "%d", ne*3) ;

  t = reshape(t, 3, ne)' ;
  
  fclose(fid) ;
