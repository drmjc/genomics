\name{miRNA.validate.miR.ID}
\alias{miRNA.validate.miR.ID}
\title{Validate miR ID's}
\usage{
  miRNA.validate.miR.ID(miRids, mirbase5p3p,
    dropped.mirs = NULL)
}
\description{
  EXPERIMENTAL: Take a vector of miR ID's, and the
  formatted miRbase table, and determine the latest miR ID
  for each supplied miRNA ALGORITHM for each miR id\cr if
  miR ID found in miRBase, then no change\cr else if miR ID
  ends in -5p or -3p, and there's only one precursor mir
  and that mir has a miR on the correct arm (5' or 3') then
  the correct miR is returned.
}
\details{
  This handles cases where miRBase uses */<> whereas query
  uses -5p/-3p else if miR ID ends in *, and there's only
  one precursor mir, and that mir has two miRs, then use
  the newest miR ID. This usually corresponds to records in
  miRbase that used to end in */<>, but now end in -5p/-3p
  else if miR ID is plain (eg miR-131) and there's only one
  precursor mir, then choose the oldest miR record. If the
  mir has only one miR then simple. if mir has 2 miRs then
  choose the oldest one. Historically mir-123 and miR-123
  would have been deposited into db when discovered, and
  later a (usually) secondard (and thus miR-123*) would
  have been discovered and deposited in db with a higher
  MIMAT accession.
}
\author{
  Mark Cowley, 2009-01-13
}

