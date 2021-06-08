
#' Make Colorscale ballanced around zero
#' Wrapper/Convenience function to gplots::colorpanel()
#' This function generates pallette of colors for a heatmap where "mid" color is centered on 0,
#' e.g. to separate the up-regulated/positive values ( "high") and down-regulated/negative values ("low")
#' @param dat_mat data matrix that will be used for ploting, should contain positive and negative values
#' @param col_res_number number of color elements in the panel/pallette, passed on to gplots::colorpanel
#' @param low optional, default "steelblue"
#' @param mid optional, default "white"
#' @param high optional, default "red"
#' @return Vector of HTML-style RGB colors
#' @export
make_balanced_colorscale = function(dat_mat , col_res_number = 10,
                                     low = "steelblue", mid = "white",  high = "red" ){
  color_balance = (max(dat_mat)  +  min(dat_mat) ) / ( max(dat_mat)  -  min(dat_mat) )
  balance_corrector = round( color_balance * col_res_number)
  if (balance_corrector < 0){
    col = gplots::colorpanel(col_res_number, low = low, mid = mid,  high = high)[ 1 : (col_res_number + balance_corrector)]
  }else{
    col = gplots::colorpanel(col_res_number, low = low, mid = mid,  high = high)[ balance_corrector : col_res_number ]
  }
  col
}
