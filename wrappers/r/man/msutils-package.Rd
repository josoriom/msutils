\name{msutils-package}
\alias{msutils}
\alias{msutils-package}
\title{msutils}
\description{
R wrappers for the \code{msutils} library for LC-MS data processing.
}
\details{
\code{msutils} allows to load mzML or ion files, extract EICs, pick peaks, and detect features across one or many samples.

A \emph{sample file} is a pointer to a parsed mzML held in memory, or to an ion file decoded on demand. The memory is freed automatically when the sample is no longer in use. Optionally call \code{\link{dispose}()} to free it manually if needed during processing.

\emph{XY data} is a pair of numeric vectors, where \code{x} is retention time and \code{y} is intensity.
}
\section{Install}{
\code{remotes::install_github("josoriom/msutils", subdir = "wrappers/r")}
}
\section{Load a file}{
\code{\link{parse_mzml}()} — load an mzML into a sample file.

\code{\link{parse_ion}()} — load an ion binary into a sample file.

\code{\link{mzml_to_ion}()} — convert a mzML to compressed ion bytes.

\code{\link{dispose}()} — free a sample file immediately.
}
\section{Convert}{
\code{\link{ion_to_json}()} — get JSON from a sample file.

\code{\link{ion_to_mzml}()} — get mzML text from a sample file.

\code{\link{ion_to_df}()} — get R lists and data frames from a sample file.
}
\section{Extract}{
\code{\link{calculate_eic}()} — get an EIC from a sample file for one m/z.

\code{\link{get_scans}()} — get scans from a sample file by RT or m/z query.

\code{\link{calculate_baseline}()} — estimate a baseline for XY data.
}
\section{Find peaks}{
\code{\link{find_peaks}()} — find all peaks in XY data.

\code{\link{get_peak}()} — find one peak near a target RT in XY data.

\code{\link{get_peaks_from_eic}()} — find peaks for many targets from a sample file.

\code{\link{get_peaks_from_chrom}()} — find peaks from stored chromatograms in a sample file.
}
\section{Find features}{
\code{\link{find_feature}()} — find targeted features by m/z and RT from a sample file.

\code{\link{find_features}()} — find all features in a sample file.

\code{\link{get_features}()} — find and align features across many sample files.
}
\keyword{package}