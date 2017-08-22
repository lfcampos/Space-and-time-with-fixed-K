#include <Rcpp.h>
using namespace Rcpp;

// CONFIRMED Jan 2017 to match Vinay smoothed PSF (psf_image_function2.R)
// up to arithematic precision

// [[Rcpp::export]]
int psf_grid_cpp(double one_point, NumericVector psf_range, NumericVector psf_grid)
{
  if ((one_point < psf_range[0]) || (one_point > psf_range[1])){
    return(-1);
  } else {
    int index = which_min(pow(one_point-psf_grid,2));
    double closest = psf_grid[index];
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
    x_indexes[i] = psf_grid_cpp(x[i],xpsf_range,xpsf_grid);
    y_indexes[i] = psf_grid_cpp(y[i],ypsf_range,ypsf_grid);
  }
  NumericVector value(m);
  double interpx1;
  double interpx2;
  double interpy1;
  double interpy2;
  double c_common;
  double c1;
  double c2;
  double c3;
  double c4;
  for (int i=0; i<m; i++){
    if((x_indexes[i]>-1) && (y_indexes[i]>-1)){
      interpx1 = xpsf_grid(x_indexes[i]);
      interpx2 = xpsf_grid(x_indexes[i]+1);
      interpy1 = ypsf_grid(y_indexes[i]);
      interpy2 = ypsf_grid(y_indexes[i]+1);
      c_common = (interpx2-interpx1)*(interpy2-interpy1);
      c1 = (interpx2-x[i])*(interpy2-y[i]);
      c2 = (x[i]-interpx1)*(interpy2-y[i]);
      c3 = (interpx2-x[i])*(y[i]-interpy1);
      c4 = (x[i]-interpx1)*(y[i]-interpy1);
      value[i] = (1/c_common)*(c1*psf_image_values(x_indexes[i],y_indexes[i])+c2*psf_image_values(x_indexes[i]+1,y_indexes[i])+c3*psf_image_values(x_indexes[i],y_indexes[i]+1)+c4*psf_image_values(x_indexes[i]+1,y_indexes[i]+1));
    }
  }
  return value;      
}












