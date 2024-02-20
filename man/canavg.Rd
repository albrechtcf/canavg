\name{hello}
\alias{hello}
\title{Hello, World!}
\usage{
hello()
}
\description{
Computes coefficients for weighted averaging of concave hemispherical mirror canopy densiometer measurements taken in each cardinal direction. This is accomplished by reading in the output of Gap Light Analyzer and calibrating mirror measurements to binned sky regions.

Arguments:
gla_detailed_report_txt = Filepath to the text file output of Gap Light Analyzer. This text output should not be modified in any way.
eye_y = The Y dimension of the point from which the densiometer is viewed. Note that the X dimension is left-right on the densiometer surface, and Y dimension is bottom-top, and the Z dimension is a line going through the middle of the engraved densiometer grid an into the ground immediately below. This defaults to -8. Units are inches.
eye_z = The Z dimension of the point from which the densiometer is viewed. Note that the X dimension is left-right on the densiometer surface, and Y dimension is bottom-top, and the Z dimension is a line going through the middle of the engraved densiometer grid an into the ground immediately below. This defaults to +10. Units are inches.
display_diagnostics = Display diagnostic graphics? Defaults to TRUE.
}
\examples{
hello()
}
