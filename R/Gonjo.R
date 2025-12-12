#' @title AMS data of the Gonjo Basin
#' @description 
#' This is the anisotropy of magnetic susceptibility (AMS) data from a 3km think section of redbeds in the Gonjo Basin in eastern Tibet that was analysed by \insertCite{li2020an;textual}{TFORGE}.
#' @references \insertAllCited{}
#' @source \doi{10.5281/zenodo.3666760}
#' @format 
#' A list with entry `datatable` containing one row per specimen and entry `matrices` containing the AMS tensor (i.e. symmetric matrix) for each specimen.
#' The `datatable` entry has 542 rows and 25 variables:
#'
#' - `Name`: character — Specimen name
#' - `really depth`: double — Depth of specimen (in meters)
#' - `Field`: double — Unsure — \insertCite{li2020an;textual}{TFORGE} say they applied a 300 A/m magnetic field
#' - `Freq.`: double — Frequency of oscillation of the applied magnetic field
#' - `Km`: double — Mean magnetic susceptibility
#' - `L`: double — Lineation (\eqn{\lambda_1 / \lambda_2})
#' - `F`: double — Foliation (\eqn{\lambda_2 / \lambda_3})
#' - `P`: double — Uncorrected degree of anisotropy (\eqn{\lambda_1 / \lambda_3} )
#' - `Pj`: double — Corrected degree of anisotropy
#' - `T`: double — Shape factor
#' - `U`: double — *Unsure*
#' - `Q`: double — *Unsure*
#' - `E`: double — *Unsure*
#' - `K1decI` and `K1incI`: doubles — In-situ direction of first eigenvector
#' - `K2decI` and `K2incI`: doubles — In-situ direction of second eigenvector
#' - `K3decI` and `K3incI`: doubles — In-situ direction of third eigenvector
#' - `K1decT` and `K1incT`: doubles — Tilt-corrected direction of first eigenvector
#' - `K2decT` and `K2incT`: doubles — Tilt-corrected of second eigenvector
#' - `K3decT` and `K3incT`: doubles — Tilt-corrected of third eigenvector
#' @details 
#' The AMS matrices were calculated using the in-situ directions by Dr. Janice Scealy.
#' 
#' The data from \doi{10.5281/zenodo.3666760} has a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/) license.
"Gonjo"

