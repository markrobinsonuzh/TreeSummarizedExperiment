context("MicrobiomeExperiment-class")

test_that("MicrobiomeExperiment-class", {
    expect_true(validObject(MicrobiomeExperiment()))
    expect_true(validObject(MicrobiomeExperiment(SimpleList())))
    sampleNames <- letters[1:4]
    cd <- DataFrame(a = letters[1:4], b = 1:4)
    rd <- DataFrame(Kingdom = "Bacteria",
                    Phylum = c("Firmicutes","Firmicutes","Firmicutes","Bacteroidetes","Euryarchaeota"),
                    Class = c("Bacilli","Bacilli","Bacilli","Saprospirae","Methanobacteria"),
                    Order = c("Bacillales","Lactobacillales","Lactobacillales","Saprospirales","Methanobacteriales"),
                    Family = c("Planococcaceae","Enterococcaceae","Enterococcaceae","Chitinophagaceae","Methanobacteriaceae"),
                    Genus = c("Staphylococcus",NA,"Melissococcus","Sediminibacterium","Methanobrevibacter"))
    counts <- matrix(sample(1:1000, nrow(rd) * nrow(cd), replace=TRUE),
                     nr = nrow(rd),
                     nc = nrow(cd))
    refSeq <- DNAStringSetList(one = DNAStringSet(c("A","A","A","A","A")),
                               two = DNAStringSet(c("A","A","A","A","A")))
    expect_error(MicrobiomeExperiment(assays = SimpleList(counts = counts),
                                      rowData = taxa,
                                      colData = cd,
                                      referenceSeq = refSeq[list(1:5,1:4)]))
    expect_error(MicrobiomeExperiment(assays = SimpleList(counts = counts),
                                      rowData = taxa,
                                      colData = cd,
                                      referenceSeq = refSeq[list(1:4,1:4)]))
    expect_error(MicrobiomeExperiment(assays = SimpleList(counts = counts),
                                      rowData = taxa,
                                      colData = cd,
                                      referenceSeq = refSeq[[1L]][1:4]))
    me <- MicrobiomeExperiment(assays = SimpleList(counts = counts),
                               rowData = rd,
                               colData = cd,
                               referenceSeq = refSeq)
    expect_error(referenceSeq(me) <- refSeq[list(1:5,1:4)])
    expect_error(referenceSeq(me) <- refSeq[list(1:4,1:4)])
    expect_error(referenceSeq(me) <- refSeq[[1L]][1:4])
    #
    se <- as(me,"SummarizedExperiment")
    expect_s4_class(as(se, "MicrobiomeExperiment"),
                  "MicrobiomeExperiment")
    expect_s4_class(as(as(se,"SingleCellExperiment"),
                       "MicrobiomeExperiment"),
                    "MicrobiomeExperiment")
    expect_s4_class(as(as(se,"TreeSummarizedExperiment"),
                       "MicrobiomeExperiment"),
                    "MicrobiomeExperiment")
    expect_s4_class(as(se, "MicrobiomeExperiment"),
                    "MicrobiomeExperiment")
    #
    referenceSeq(me) <- refSeq
    expect_s4_class(referenceSeq(me),"DNAStringSetList")
    referenceSeq(me) <- as.list(refSeq)
    expect_s4_class(referenceSeq(me),"DNAStringSetList")
    referenceSeq(me) <- refSeq[[1L]]
    expect_s4_class(referenceSeq(me),"DNAStringSet")
    referenceSeq(me) <- as.character(refSeq[[1L]])
    expect_s4_class(referenceSeq(me),"DNAStringSet")
})
