function [x,y,z,d]=readbbox(file)

  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%f") ;
  nx = length(dat)/7 ;
  dat = reshape(dat, 7, nx)' ;

  d = dat(:, 7) ;

  x = dat(:, [1 2 2 1 1]) ;
  y = dat(:, [3 3 4 4 3]) ;
  z = dat(:, [5 5 6 6 5]) ;
  
  fclose(fid) ;
