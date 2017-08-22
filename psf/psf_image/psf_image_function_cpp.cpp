#include <Rcpp.h>
using namespace Rcpp;

// CONFIRMED Jan 2017 to match Vinay grid PSF (psf_image_function.R)

// [[Rcpp::export]]
int psf_gridx_cpp(double one_point, NumericVector xpsf_range, NumericVector xpsf_grid)
{
  if ((one_point < xpsf_range[0]) || (one_point > xpsf_range[1])){
    return(-1);
  } else {
    int index = which_min(pow(one_point-xpsf_grid,2));
    double closest = xpsf_grid[index];
    if (one_point < closest){
      return(index-1);
    } else {
      return(index);
    }
  }  
}

// [[Rcpp::export]]
int psf_gridy_cpp(double one_point, NumericVector ypsf_range, NumericVector ypsf_grid)
{
  if ((one_point < ypsf_range[0]) || (one_point > ypsf_range[1])){
    return(-1);
  } else {
    int index = which_min(pow(one_point-ypsf_grid,2));
    double closest = ypsf_grid[index];
    if (one_point < closest){
      return(index-1);
    } else {
      return(index);
    }
  }  
}

// [[Rcpp::export]]
NumericVector psf_image_cpp(NumericVector pts, NumericVector center, NumericMatrix psf_image_values, NumericVector xpsf_range, NumericVector xpsf_grid, NumericVector ypsf_range, NumericVector ypsf_grid)
{
  int m = pts.size()/2;
  double x[m]; double y[m];
  for (int i=0; i<m; i++){
    x[i] = pts[i] - center[0];
    y[i] = pts[i+m] - center[1];
  }
  NumericVector x_indexes(m);
  NumericVector y_indexes(m);
  for (int i=0; i<m; i++){
    x_indexes[i] = psf_gridx_cpp(x[i],xpsf_range,xpsf_grid);
    y_indexes[i] = psf_gridy_cpp(y[i],ypsf_range,ypsf_grid);
  }
  NumericVector value(m);
  for (int i=0; i<m; i++){
    if((x_indexes(i)>-1) && (y_indexes(i)>-1)){
      value(i) = psf_image_values(x_indexes(i),y_indexes(i));
    }
  }
  return value;      
}












