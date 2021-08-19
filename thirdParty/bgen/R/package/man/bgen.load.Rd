% File src/library/base/man/load.Rd
\name{bgen.load}
\alias{bgen.load}
\alias{bgen}
\alias{rbgen}
\title{Load genotype data from an indexed BGEN file}
\author{Gavin Band}
\description{
Loads genotypes and associated metadata from specified regions, or specified variants, in a BGEN file.
The BGEN file must have been previously indexed using bgenix.
}
\usage{
bgen.load(
  filename,
  ranges = NULL,
  rsids = NULL,
  max_entries_per_sample = 3,
  samples = NULL,
  index.filename = sprintf( "\%s.bgi", filename )
)
}
\arguments{
  \item{filename}{The name of a BGEN file to load data from.  The corresponding bgenix index file (file.bgi) must also exist.}
  \item{ranges}{A dataframe with 'chromosome', 'start', and 'end' columns, holding genomic regions for which to load data.}
  \item{rsids}{A character vector holding the IDs of variants for which to load data.}
  \item{max_entries_per_sample}{An integer specifying the maximum number of probabilities expected per variant per sample.  This is used to set the third dimension of the data matrix returned.  For biallelic variants and diploid samples, there are three possible genotypes at each variant so the default value will work. For more complex situations you may need to increase this to allocate sufficient space.}
  \item{samples}{A character vector specifying the IDs of samples to load data for.}
  \item{index.filename}{The name of the index file, if different from the default.}
}
\details{
First, the specified BGEN file and its corresponding index file are opened.  Then the index file
is used to locate the data for the specified ranges and rsids.  Then the genotype data is loaded.
}
\value{
A list with several members:
 \item{variants}{a data frame containing the genomic position, IDs, and alleles of loaded variants.  These are ordered as specified by the index, which will usually be order of genomic position.}
 \item{samples}{a vector of identifiers of the samples for which data has been loaded.}
 \item{ploidy}{a matrix of ploidy values, with one row per variant and one column per sample.}
 \item{phased}{a vector of logical values indicating whether data for the variants is phased or not}
 \item{data}{an array of genotype probability values, indexed by variants, samples, and by genotype.}
}
\seealso{
}
\examples{
D = bgen.load( "example/example.16bits.bgen", rsids = c( "RSID_101", "RSID_40" ))
print( D$data[,'sample_001',] )
print( D$data[ 'RSID_40',, ] )
bgen.load( "example/example.16bits.bgen", rsids = c( "RSID_101", "RSID_40" ), samples = c( "sample_001", "sample_101"))
bgen.load( "example/example.16bits.bgen", ranges = data.frame( chromosome = '01', start = 0, end = 5000 ))
D = bgen.load( "example/complex.bgen", ranges = data.frame( chromosome = '01', start = 0, end = 1000 ), max_entries_per_sample = 40 )
D$data[4,,]
}
\keyword{bgen}
