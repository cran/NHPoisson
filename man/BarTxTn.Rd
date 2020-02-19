\name{BarTxTn}
\alias{BarTxTn}
\docType{data}
\title{Barcelona temperature data}

\description{
Barcelona daily temperature series during the summer months (May, June, July, 
August and September) from  1951 to 2004.}

\usage{data(BarTxTn)}


\details{Variables

dia: Postion of the day in the year, from 121 (1st of May) to 253 (30th of September).

mes: Month of the year, from 5 to 9.
 
ano: Year, from 1951 to 2004.

diames: Position of the day in the month, from 1 to 30 or 31.

Tx: Daily maximum temperature.

Tn: Daily minimum temperature.

Txm31: Local maximum temperature signal. Lowess of Tx with a centered window of 31 days.

Txm15: Local maximum temperature signal. Lowess of Tx with a centered window of 15 days.

Tnm31: Local minimum temperature signal. Lowess of Tn with a centered window of 31 days.

Tnm15: Local minimum temperature signal. Lowess of Tn with a centered window of 15 days.

TTx: Long term maximum temperature signal. Lowess of Tx with a centered 40\% window.

TTn: Long term minimum temperature signal. Lowess of Tn with a centered 40\% window.

}


\references{


Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

}


\examples{
data(BarTxTn)
}

