\name{msutils-package}
\alias{msutils}
\alias{msutils-package}
\title{msutils}
\description{
Simple R wrappers for the \code{msutils} Rust library.
}
\details{
\code{msutils} lets you load mzML or binary files, get EICs, find peaks,
find features, and align features across many samples.
}
\section{Install}{
\code{remotes::install_github("josoriom/msutils", subdir = "wrappers/r")}
}
\section{Main functions}{
\code{\link{parse_mzml}()} reads mzML bytes and returns a file handle.

\code{\link{parse_bin}()} reads a binary msutils file and returns a file handle.

\code{\link{mzml_to_bin}()} converts a loaded file into binary bytes.

\code{\link{calculate_eic}()} gets an EIC for one target m/z.

\code{\link{calculate_baseline}()} estimates a baseline for one signal.

\code{\link{find_peaks}()} finds peaks in one trace.

\code{\link{get_peak}()} gets one target peak near one RT.

\code{\link{get_peaks_from_chrom}()} finds many target peaks from chromatograms.

\code{\link{get_peaks_from_eic}()} finds many target peaks from EICs.

\code{\link{find_feature}()} finds one or more target features.

\code{\link{find_features}()} finds untargeted features in one file.

\code{\link{get_scans}()} gets scans by RT range, closest RT, m/z range, or closest m/z.

\code{\link{get_features}()} aligns features across all files in one folder.

\code{\link{bin_to_json}()}, \code{\link{bin_to_mzml}()}, and
\code{\link{bin_to_df}()} convert a loaded file into other R outputs.

\code{\link{dispose}()} frees a file handle when you want to release memory now.
}
\keyword{package}
