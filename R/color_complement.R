
#' Create color complement by rotating the color hue
#'
#' Create color complement by rotating the color hue
#'
#' This function rotates the color hue to create a complementary
#' color for each `color` input. It differs from standard methods
#' by using warped color hue by default (`useWarpHue=TRUE`), which
#' uses a red-yellow-blue color wheel instead of R default
#' red-green-blue. It also imposes a minimum chroma, which
#' ensures the output color is reasonably high in color
#' saturation.
#'
#' @family jam utility functions
#'
#' @param color `character` vector of R compatible colors.
#' @param Hflip `numeric` value in degrees (from 0 to 360) added
#'    to the color hue to produce the final color hue. Typically
#'    180 degrees will select the color opposite the input color
#'    on a virtual color wheel. Note that `warpHue=TRUE` also
#'    enables a customized color wheel.
#' @param Cfloor `numeric` value used to limit output chroma `C`
#'    values to this minimum value, to ensure a minimum color saturation.
#' @param Lrange `numeric` vector with the allowed range of output
#'    luminance `L` values. When supplied, output values are
#'    simply forced to this range with no other scaling of intermediate
#'    values.
#' @param useWarpHue `logical` indicating whether to use the warp
#'    hue functions `colorjam::h2hw()` and `colorjam::hw2h()` which
#'    effectively change the color wheel from red-green-blue to
#'    red-yellow-blue.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' n <- 5;
#' rc <- colorjam::rainbowJam(n);
#' rc_comp <- color_complement(rc);
#' rc_comp2 <- color_complement(rc, useWarpHue=FALSE);
#' jamba::showColors(list(rainbowJam=rc,
#'    `complement\n(preset="dichromat")`=rc_comp,
#'    `complement\n(useWarpHue=FALSE)`=rc_comp2));
#'
#' rc <- colorjam::rainbowJam(n, preset="ryb");
#' rc_comp <- color_complement(rc, preset="ryb");
#' rc_comp2 <- color_complement(rc, useWarpHue=FALSE);
#' jamba::showColors(list(`rainbowJam\n(preset="ryb")`=rc,
#'    `complement\n(preset="ryb")`=rc_comp,
#'    `complement\n(useWarpHue=FALSE)`=rc_comp2));
#'
#' ## divergent color gradients through white
#' ## hint: use higher lens value to make middle colors more intense
#' rc <- colorjam::rainbowJam(n);
#' rc_comp <- color_complement(rc);
#' rc_ramps <- lapply(jamba::nameVector(seq_along(rc)), function(i){
#'    j <- jamba::getColorRamp(c(rc[i], "white", rc_comp[i]),
#'       n=25,
#'       lens=0,
#'       divergent=TRUE);
#'    names(j) <- "";
#'    names(j)[1] <- "original colors";
#'    names(j)[25] <- "color complements";
#'    j;
#' });
#' jamba::showColors(rc_ramps, groupCellnotes=TRUE, groupByColors=FALSE);
#'
#' ## divergent color gradients through black
#' ## hint: use higher lens value to make middle colors more intense
#' rc_ramps2 <- lapply(jamba::nameVector(seq_along(rc)), function(i){
#'    j <- jamba::getColorRamp(c(rc[i], "black", rc_comp[i]),
#'       n=25,
#'       lens=1,
#'       divergent=TRUE);
#'    names(j) <- "";
#'    names(j)[1] <- "original colors";
#'    names(j)[25] <- "color complements";
#'    j;
#' });
#' jamba::showColors(rc_ramps2, groupCellnotes=TRUE, groupByColors=FALSE);
#'
#' @export
color_complement <- function
(color,
 Hflip=180,
 Cfloor=NULL,
 Crange=c(80, 100),
 Lrange=c(50, 85),
 preset="dichromat",
 useWarpHue=TRUE,
 ...)
{
   if (length(jamba::rmNA(Crange)) == 0) {
      Crange <- c(60, 100);
   } else {
      Crange <- range(jamba::rmNA(Crange));
   }
   if (length(Cfloor) == 0) {
      Cfloor <- min(Crange, na.rm=TRUE);
   }

   # convert input to HCL
   hcl <- jamba::col2hcl(color);
   H <- hcl["H",];
   # optionally adjust color wheel
   if (TRUE %in% useWarpHue) {
      H <- colorjam::h2hw(h=H,
         preset=preset);
   }
   newH <- (Hflip + H) %% 360;
   # optionally revert the adjusted color wheel
   if (TRUE %in% useWarpHue) {
      newH <- colorjam::hw2h(h=newH,
         preset=preset);
   }
   hcl["H",] <- newH;

   # apply Crange
   if (any(hcl["C",] < Cfloor | hcl["C",] > max(Crange, na.rm=TRUE))) {
      hcl["C",] <- jamba::noiseFloor(
         hcl["C",],
         ceiling=max(Crange, na.rm=TRUE),
         minimum=Cfloor);
   }
   # apply Lrange
   if (length(Lrange) > 0) {
      hcl["L",] <- jamba::noiseFloor(
         hcl["L",],
         minimum=min(Lrange, na.rm=TRUE),
         ceiling=max(Lrange, na.rm=TRUE));
   }
   # convert back to hex
   color2 <- jamba::hcl2col(hcl);
   return(color2);
}
