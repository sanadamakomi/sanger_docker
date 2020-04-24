## time
get_time_human <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%OS")
}

## fig height
figheight <- function(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE) {
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(aveposition)] #####
    basecalls2 <- basecalls2[1:length(aveposition)] ######
    if(showtrim == FALSE) {
        if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
        else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
        if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
        else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
        aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
    }
    indexes <- 1:length(basecalls1)
    trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all false if not trimmed
    if (!is.null(trim3)) {
        traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, nrow(traces))),]
    }
    if (!is.null(trim5)) {
        offset <- max(c(1, aveposition[1] - 10))
        traces <- traces[offset:nrow(traces),]
        aveposition <- aveposition - (offset-1)
    }
  
    valuesperbase <- nrow(traces)/length(basecalls1)
    tracewidth <- width*valuesperbase
    breaks <- seq(1,nrow(traces), by=tracewidth) 
  
    numplots <- length(breaks)
    return(numplots*pixelsperrow)
}

## clean string
cleanString <- function(string) {
    string <- gsub(" +", "", string)
    string <- gsub("\n?", "", string)
    return(string)
}

## format position
formatPos <- function(string) {
    string <- cleanString(string)
    if (string == "") return(NULL)
    pos <- gsub("^CHR", "", toupper(string), ignore.case = TRUE)
    pos1 <- str_match(pos, "^(.+):([0-9]+)-([0-9]+)$")
    pos2 <- str_match(pos, "^(.+):([0-9]+)$")
    if (nrow(na.omit(pos1)) == 0 & nrow(na.omit(pos2)) == 0) {
        return(NULL)
    }
    if (nrow(na.omit(pos1)) == 1) {
        s <- as.numeric(pos1[1,3])
        e <- as.numeric(pos1[1,4])
        if (s <= e) {
            return(list(chr=pos1[1,2], check.start = s, check.end = e))
        } else {
            return(list(chr=pos1[1,2], check.start = e, check.end = s))
        }        
    }
    if (nrow(na.omit(pos2)) == 1) {
    return(list(chr=pos2[1,2], check.start = as.numeric(pos2[1,3]), check.end = as.numeric(pos2[1,3])))
    }
}

## creat seqinfo object
getSeqinfo <- function() {
    seq.vec <- c("1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663, "8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540, "15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520, "21"=48129895, "22"=51304566, "X"=155270560, "Y"=59373566, "MT"=16569)
    return(Seqinfo(seqnames = names(seq.vec), seqlengths = seq.vec, isCircular=rep(NA, length(seq.vec)), genome=rep(NA, length(seq.vec))))
}

## read Blast result
readBlast <- function(text) {
    indat <- try(read.table(text = text, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(indat, 'try-error')) {
        write("Warning: missing blast", stdout())
        return(NULL)
    }
    if (nrow(indat) > 0) {
        colnames(indat) <- c("query.id", "subject.id", "identity", "alignment.length", "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "qseq", "sseq", "evalue", "bit.score", "subject.tax.ids", "subject.titles")
        indat[,c("query.id", "subject.id")] <- apply(indat[,c("query.id", "subject.id"), drop=FALSE], 2, function(x) {as.character(as.vector(x))})
        indat[,c("identity", "q.start", "q.end", "s.start", "s.end")] <- apply(indat[,c("identity", "q.start", "q.end", "s.start", "s.end"), drop=FALSE], 2, function(x) {as.character(as.vector(x))})
        indat[,"subject.id"] <- gsub(" +", "", as.character(as.vector(indat[,"subject.id"])))
        idx <- which(as.character(as.vector(indat[,"subject.id"])) %in% c(as.character(seq(1,22,1)), "X", "Y", "MT"))
        if (length(idx) > 0) return(indat[idx,, drop = FALSE])
    }
    write("Warning: missing blast", stdout())
    return(NULL)
}

## seq to blast result
doBlast <- function(seq, toolDir, dbPrefix) {
    ## do blast
    seq <- as.character(seq)
    cmd <- paste0("echo -e \">sequence\n", seq, "\" | ", toolDir, "/blastn -db ", dbPrefix, " -max_target_seqs 3 -max_hsps 3 -outfmt ", shQuote("7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore staxids salltitles", type = "cmd"))
    readBlast(system(cmd, intern = TRUE))
}

## blast result to gr
blastToGR <- function(x) {
    seq.info <- getSeqinfo()
    chr <- as.character(as.vector(x[,"subject.id"]))
    ss <- as.numeric(x[, "s.start"])
    se <- as.numeric(x[, "s.end"])
    idx <- which(ss > se)
    ss[idx] <- as.numeric(x[idx, "s.end"])
    se[idx] <- as.numeric(x[idx, "s.start"])
    GRanges(Rle(chr), IRanges(start = ss, end = se), seqinfo = seq.info)    
}

## align blast results' seq one by one between query and subject
alignMatrix <- function(x, sanger.length) {
    chr <- as.character(as.vector(x[1,"subject.id"]))
    q.start <- as.numeric(x[1, "q.start"])
    q.end <- as.numeric(x[1, "q.end"])
    s.start <- as.numeric(x[1, "s.start"])
    s.end <- as.numeric(x[1, "s.end"])
    strand <- ifelse(s.start <= s.end, 1, -1)
    q <- as.character(x[1, "qseq"])
    s <- as.character(x[1, "sseq"])
    qs <- unlist(strsplit(q, ""))
    ss <- unlist(strsplit(s, ""))
    ma <- matrix(data = c(qs, ss, rep(NA, length(qs)*2)), ncol = 4, byrow = FALSE)
    q.counter <- q.start - 1
    s.counter <- s.start - strand
    for (i in 1:nrow(ma)) {
        if (! ma[i,1] == "-") {
            q.counter = q.counter + 1
            ma[i,3] = q.counter
        } else {
            ma[i,3] = q.counter
        }
        if (! ma[i,2] == "-") {
            s.counter = s.counter + strand
            ma[i,4] = s.counter
        } else {
            ma[i,4] = s.counter
        }
    }
    trim5 <- as.numeric(ma[,3]) - 1
    trim3 <- sanger.length - as.numeric(ma[,3])
    ma <- cbind(ma, rep(strand, nrow(ma)), rep(chr, nrow(ma)), 1:nrow(ma), trim5, trim3, rep(q.start, nrow(ma)), rep(q.end, nrow(ma)), ma[,1], ma[,2], rep(q, nrow(ma)), rep(s,nrow(ma)))
    if (strand < 0) {
        ma[,1] <- chartr("ATCG", "TAGC", ma[,1])
        ma[,2] <- chartr("ATCG", "TAGC", ma[,2])
    }    
    colnames(ma) <- c("alt", "ref", "query.loc", "subject.loc", "strand", "chr", "query.loc.prime", "trim5", "trim3", "query.start", "query.end", "alt.prime", "ref.prime", "query.seq.prime", "subject.seq.prime")
    ma
}

# marge the same ref pos variants
mergeVariantMa <- function(x) {
    qt <- table(x$query.loc)
    st <- table(x$subject.loc)
    qover1 <- names(qt)[which(qt > 1)]
    sover1 <- names(st)[which(st > 1)]
    idxover1 <- c()
    out.q <- NULL
    out.s <- NULL
    out.other <- NULL
    if (length(qover1) > 0) {       
        qidxover1 <- which(x$query.loc %in% qover1)
        idxover1 <- c(idxover1, qidxover1)
        out.q <- lapply(qover1, function(p) {
            idx <- which(x$query.loc == p)
            alt <- paste(x[idx, "alt"], collapse="")
            ref <- paste(x[idx, "ref"], collapse="")
            if (x[idx[1],"strand"] < 0) {
                alt <- reverse(alt)
                ref <- reverse(ref)
            }
            alt <- gsub("-+", "", alt)
            ref <- gsub("-+", "", ref)
            if (alt == "") alt <- "-"
            if (ref == "") ref <- "-"            
            df <- data.frame(
                alt=alt,
                ref=ref,
                query.loc.start=min(x[idx, "query.loc"]),
                query.loc.end=max(x[idx, "query.loc"]),
                subject.loc.start=min(x[idx, "subject.loc"]),
                subject.loc.end=max(x[idx, "subject.loc"]),
                strand=x[idx[1], "strand"],
                chr=x[idx[1], "chr"],
                query.loc.prime.start=min(x[idx, "query.loc.prime"]),
                query.loc.prime.end=max(x[idx, "query.loc.prime"]),
                trim5=min(x[idx, "trim5"]),
                trim3=min(x[idx, "trim3"]),
                query.start=x[idx[1],"query.start"],
                query.end=x[idx[1], "query.end"],
                alt.prime=paste(x[idx, "alt.prime"], collapse=""),
                ref.prime=paste(x[idx, "ref.prime"], collapse=""),
                query.seq.prime=x[idx[1], "query.seq.prime"],
                subject.seq.prime=x[idx[1], "subject.seq.prime"]
            )
       })
        out.q <- do.call("rbind", out.q)
    }
    
    if (length(sover1) > 0) {
        sidxover1 <- which(x$subject.loc %in% sover1)
        idxover1 <- c(idxover1, sidxover1)
        out.s <- lapply(sover1, function(p) {
            idx <- which(x$subject.loc == p)
            alt <- paste(x[idx, "alt"], collapse="")
            ref <- paste(x[idx, "ref"], collapse="")
            if (x[idx[1],"strand"] < 0) {
                alt <- reverse(alt)
                ref <- reverse(ref)
            }
            alt <- gsub("-+", "", alt)
            ref <- gsub("-+", "", ref)
            if (alt == "") alt <- "-"
            if (ref == "") ref <- "-"            
            df <- data.frame(
                alt=alt,
                ref=ref,
                query.loc.start=min(x[idx, "query.loc"]),
                query.loc.end=max(x[idx, "query.loc"]),
                subject.loc.start=min(x[idx, "subject.loc"]),
                subject.loc.end=max(x[idx, "subject.loc"]),
                strand=x[idx[1], "strand"],
                chr=x[idx[1], "chr"],
                query.loc.prime.start=min(x[idx, "query.loc.prime"]),
                query.loc.prime.end=max(x[idx, "query.loc.prime"]),
                trim5=min(x[idx, "trim5"]),
                trim3=min(x[idx, "trim3"]),
                query.start=x[idx[1],"query.start"],
                query.end=x[idx[1], "query.end"],
                alt.prime=paste(x[idx, "alt.prime"], collapse=""),
                ref.prime=paste(x[idx, "ref.prime"], collapse=""),
                query.seq.prime=x[idx[1], "query.seq.prime"],
                subject.seq.prime=x[idx[1], "subject.seq.prime"]
            )
       })
        out.s <- do.call("rbind", out.s)
    }
    
    idxequal1 <- setdiff(1:nrow(x), idxover1)
    if (length(idxequal1) > 0 ) {
        out.other <- data.frame(
                alt=x[idxequal1, "alt"],
                ref=x[idxequal1, "ref"],
                query.loc.start=x[idxequal1, "query.loc"],
                query.loc.end=x[idxequal1, "query.loc"],
                subject.loc.start=x[idxequal1, "subject.loc"],
                subject.loc.end=x[idxequal1, "subject.loc"],
                strand=x[idxequal1, "strand"],
                chr=x[idxequal1, "chr"],
                query.loc.prime.start=x[idxequal1, "query.loc.prime"],
                query.loc.prime.end=x[idxequal1, "query.loc.prime"],
                trim5=x[idxequal1, "trim5"],
                trim3=x[idxequal1, "trim3"],
                query.start=x[idxequal1,"query.start"],
                query.end=x[idxequal1, "query.end"],
                alt.prime=x[idxequal1, "alt.prime"],
                ref.prime=x[idxequal1, "ref.prime"],
                query.seq.prime=x[idxequal1, "query.seq.prime"],
                subject.seq.prime=x[idxequal1, "subject.seq.prime"]
            )        
    }
    result.lst <- list(out.q, out.s, out.other)
    result.lst <- result.lst[which(! sapply(result.lst, is.null))]
    result.lst <- do.call("rbind", result.lst)
    return(result.lst)
}

## call variant from align matrix
variantMatrix <- function(x, sanger.length, trim.bin = 0) {
    ma <- alignMatrix(x, sanger.length)
    trim.start <- trim.bin + 1
    trim.end <- sanger.length - trim.bin
    idx <- which(! ma[,1] == ma[,2])
    if (length(idx) > 0) {
        ma.variant <- ma[idx,,drop = FALSE]
        q.variant.pos <- as.numeric(ma.variant[,3])
        idx1 <- which(as.numeric(ma.variant[,3]) %in% seq(trim.start,trim.end, 1))        
        if (length(idx1) > 0) {
            ma.variant.trim <- ma.variant[idx1,,drop = FALSE]
            ma.variant.trim <- as.data.frame(ma.variant.trim)  
            ma.variant.trim[,c("query.loc", "subject.loc", "strand", "query.loc.prime", "trim5", "trim3", "query.start", "query.end")] <- apply(ma.variant.trim[,c("query.loc", "subject.loc", "strand", "query.loc.prime", "trim5", "trim3", "query.start", "query.end"),drop=FALSE], 2, function(x) {as.numeric(as.vector(x))})            
            ma.variant.trim[,c("alt", "ref", "chr", "alt.prime", "ref.prime", "query.seq.prime", "subject.seq.prime")] <- apply(ma.variant.trim[,c("alt", "ref", "chr", "alt.prime", "ref.prime", "query.seq.prime", "subject.seq.prime"),drop=FALSE], 2, function(x) {as.character(as.vector(x))})            
            ma.variant.trim.merge <- mergeVariantMa(ma.variant.trim)
            return(ma.variant.trim.merge)
        }
    }
    return(NULL)
}

## call single
callSingleVariant <- function(var.dat, seqtype, sanger.length, pos.lst, trim.bin = 0) {
    check.region <- NULL
    if (! is.null(pos.lst)) {
        check.region <- GRanges(Rle(as.character(pos.lst$chr)), IRanges(start = as.numeric(pos.lst$check.start), end = as.numeric(pos.lst$check.end)), seqinfo = getSeqinfo())
    }
    var.dat.gr <- blastToGR(var.dat)    
    tokens <- 1:nrow(var.dat)
    is.match.pos <- "."
    if (! is.null(check.region)) {
        check.hit <- findOverlaps(var.dat.gr, check.region)
        tokens <- intersect(1:nrow(var.dat), queryHits(check.hit))
        is.match.pos <- TRUE
        if (length(tokens) == 0) {
            write("Warning: No hit result in target region. Use all blast result to call variants.", stdout())
            tokens <- 1:nrow(var.dat)
            is.match.pos <- FALSE
        }
    }
    result.lst <- lapply(tokens, function(i) {
        var.ma <- variantMatrix(var.dat[i,], sanger.length = sanger.length, trim.bin = trim.bin)
        if (is.null(var.ma)) {
            return(NULL)
        }
        var.ma <- as.data.frame(var.ma)
        out <- data.frame(
            chr=var.ma$chr,
            start=var.ma$subject.loc.start,
            end=var.ma$subject.loc.end,
            ref=var.ma$ref,
            alt=var.ma$alt,
            gt=rep("HET", nrow(var.ma)),
            sanger=rep(seqtype, nrow(var.ma)),
            rank=rep(i, nrow(var.ma)),
            is.match.pos=rep(is.match.pos, nrow(var.ma)),
            strand=as.numeric(var.ma$strand),
            trim5=as.numeric(var.ma$trim5),
            trim3=as.numeric(var.ma$trim3),
            query.loc.start=as.numeric(var.ma$query.loc.start),
            query.loc.end=as.numeric(var.ma$query.loc.end),
            query.loc.prime.start=as.numeric(var.ma$query.loc.prime.start),
            query.loc.prime.end=as.numeric(var.ma$query.loc.prime.end),  
            query.start=as.numeric(var.ma$query.start),
            query.end=as.numeric(var.ma$query.end),
            ref.prime=var.ma$ref.prime,
            alt.prime=var.ma$alt.prime,
            query.seq.prime=var.ma$query.seq.prime,
            subject.seq.prime=var.ma$subject.seq.prime             
        )
        return(out)
    })
    result.lst <- result.lst[which(! sapply(result.lst, is.null))]
    if (length(result.lst) > 0 ) {
        result.lst <- do.call("rbind", result.lst)
    } else {
        result.lst <- NULL
    }
    return(result.lst)
}

## call variant
callVariant <- function(ps.dat, ss.dat, sanger.length, pos.lst, trim.bin = 0) {
    check.region <- NULL
    if (! is.null(pos.lst)) {
        check.region <- GRanges(Rle(as.character(pos.lst$chr)), IRanges(start = as.numeric(pos.lst$check.start), end = as.numeric(pos.lst$check.end)), seqinfo = getSeqinfo())
    }

    ps.dat.gr <- blastToGR(ps.dat)
    ss.dat.gr <- blastToGR(ss.dat)
    hit <- findOverlaps(ps.dat.gr, ss.dat.gr)
    tokens <- 1:nrow(ps.dat)
    is.match.pos <- "."
    if (! is.null(check.region)) {
        check.hit <- findOverlaps(ps.dat.gr, check.region)
        tokens <- intersect(1:nrow(ps.dat), queryHits(check.hit))
        is.match.pos <- TRUE
        if (length(tokens) == 0) {
            write("Warning: No hit result in target region. Use all blast result to call variants.", stdout())
            tokens <- 1:nrow(ps.dat)
            is.match.pos <- FALSE
        }
    }
    result.lst <- lapply(tokens, function(i) {
        if (! i %in% queryHits(hit)) {
            ps.var.ma <- variantMatrix(ps.dat[i,], sanger.length = sanger.length, trim.bin = trim.bin)
            if (is.null(ps.var.ma)) {
                return(NULL)
            }
            ps.var.ma <- as.data.frame(ps.var.ma)
            out <- data.frame(
                chr=ps.var.ma$chr,
                start=ps.var.ma$subject.loc.start,
                end=ps.var.ma$subject.loc.end,
                ref=ps.var.ma$ref,
                alt=ps.var.ma$alt,
                gt=rep("HOM", nrow(ps.var.ma)),
                sanger=rep("primary", nrow(ps.var.ma)),
                rank=rep(i, nrow(ps.var.ma)),
                is.match.pos=rep(is.match.pos, nrow(ps.var.ma)),
                strand=as.numeric(ps.var.ma$strand),
                trim5=as.numeric(ps.var.ma$trim5),
                trim3=as.numeric(ps.var.ma$trim3),
                query.loc.start=as.numeric(ps.var.ma$query.loc.start),
                query.loc.end=as.numeric(ps.var.ma$query.loc.end),
                query.loc.prime.start=as.numeric(ps.var.ma$query.loc.prime.start),
                query.loc.prime.end=as.numeric(ps.var.ma$query.loc.prime.end),  
                query.start=as.numeric(ps.var.ma$query.start),
                query.end=as.numeric(ps.var.ma$query.end),
                ref.prime=ps.var.ma$ref.prime,
                alt.prime=ps.var.ma$alt.prime,
                query.seq.prime=ps.var.ma$query.seq.prime,
                subject.seq.prime=ps.var.ma$subject.seq.prime             
            )
            return(out)
        } else {
            j <- min(subjectHits(hit[which(queryHits(hit) == i)]))
            ps.var.ma <- variantMatrix(ps.dat[i,], sanger.length = sanger.length, trim.bin = trim.bin)
            ss.var.ma <- variantMatrix(ss.dat[j,], sanger.length = sanger.length, trim.bin = trim.bin)
            if (is.null(ps.var.ma) & is.null(ss.var.ma)) {
                return(NULL)
            } else if (is.null(ps.var.ma) | is.null(ss.var.ma)) {
                if (!is.null(ps.var.ma)) {
                    var.ma <- as.data.frame(ps.var.ma)
                    sanger <- "primary"
                } else {
                    var.ma <- as.data.frame(ss.var.ma)
                    sanger <- "secondary"
                }
                out <- data.frame(
                    chr=var.ma$chr,
                    start=var.ma$subject.loc.start,
                    end=var.ma$subject.loc.end,
                    ref=var.ma$ref,
                    alt=var.ma$alt,
                    gt=rep("HET", nrow(var.ma)),
                    sanger=rep(sanger, nrow(var.ma)),
                    rank=rep(i, nrow(var.ma)),
                    is.match.pos=rep(is.match.pos, nrow(var.ma)),
                    strand=as.numeric(var.ma$strand),
                    trim5=as.numeric(var.ma$trim5),
                    trim3=as.numeric(var.ma$trim3),
                    query.loc.start=as.numeric(var.ma$query.loc.start),
                    query.loc.end=as.numeric(var.ma$query.loc.end),
                    query.loc.prime.start=as.numeric(var.ma$query.loc.prime.start),
                    query.loc.prime.end=as.numeric(var.ma$query.loc.prime.end),  
                    query.start=as.numeric(var.ma$query.start),
                    query.end=as.numeric(var.ma$query.end),
                    ref.prime=var.ma$ref.prime,
                    alt.prime=var.ma$alt.prime,
                    query.seq.prime=var.ma$query.seq.prime,
                    subject.seq.prime=var.ma$subject.seq.prime                    
                )
                return(out)
            } else {
                ps.position <- paste0(ps.var.ma$chr, ":", ps.var.ma$subject.loc.start, "-", ps.var.ma$subject.loc.end)
                ss.position <- paste0(ss.var.ma$chr, ":", ss.var.ma$subject.loc.start, "-", ss.var.ma$subject.loc.end)
                var.pos.both <- intersect(ps.position, ss.position)
                var.pos.ps <- setdiff(ps.position, ss.position)
                var.pos.ss <- setdiff(ss.position, ps.position)
                out.both <- NULL
                out.ps <- NULL
                out.ss <- NULL
                if (length(var.pos.both) > 0) {
                    out.both <- lapply(var.pos.both, function(p) {             
                        pidx.p <- which(ps.position == p)
                        pidx.s <- which(ss.position == p)
                        ref <- as.character(as.vector(ps.var.ma[pidx.p, "ref"]))
                        alt.p <- ps.var.ma[pidx.p, "alt"]
                        alt.s <- ss.var.ma[pidx.s, "alt"]
                        alt.p <- paste(alt.p, collapse = "")
                        alt.s <- paste(alt.s, collapse = "")
                        if (alt.p == alt.s) {
                            alt <- alt.p
                            gt <- "HOM"
                        } else {
                            alt <- paste(alt.p, alt.s, sep = ",")
                            gt <- "HET"
                        }                        
                        df <- data.frame(
                            chr=as.character(as.vector(ps.var.ma[pidx.p, "chr"])),
                            start=ps.var.ma[pidx.p, "subject.loc.start"],
                            end=ps.var.ma[pidx.p, "subject.loc.end"],
                            ref=ref,
                            alt=alt,
                            gt=gt,
                            sanger="both",
                            rank=i,
                            is.match.pos=is.match.pos,
                            strand=ps.var.ma[pidx.p, "strand"],
                            trim5=ps.var.ma[pidx.p, "trim5"],
                            trim3=ps.var.ma[pidx.p, "trim3"],
                            query.loc.start=ps.var.ma[pidx.p, "query.loc.start"],
                            query.loc.end=ps.var.ma[pidx.p, "query.loc.end"],
                            query.loc.prime.start=ps.var.ma[pidx.p, "query.loc.prime.start"],
                            query.loc.prime.end=ps.var.ma[pidx.p, "query.loc.prime.end"],
                            query.start=ps.var.ma[pidx.p, "query.start"],
                            query.end=ps.var.ma[pidx.p, "query.end"],
                            ref.prime=ps.var.ma[pidx.p, "ref.prime"],
                            alt.prime=ps.var.ma[pidx.p, "alt.prime"],
                            query.seq.prime=ps.var.ma[pidx.p, "query.seq.prime"],
                            subject.seq.prime=ps.var.ma[pidx.p, "subject.seq.prime"]
                        )
                        return(df)
                    })
                    out.both <- do.call("rbind", out.both)                
                }
                if (length(var.pos.ps) > 0) {
                    pidx.p <- which(ps.position %in% var.pos.ps)
                    out.ps <- data.frame(
                        chr=as.character(as.vector(ps.var.ma[pidx.p, "chr"])),
                        start=ps.var.ma[pidx.p, "subject.loc.start"],
                        end=ps.var.ma[pidx.p, "subject.loc.end"],
                        ref=ps.var.ma[pidx.p, "ref"],
                        alt=ps.var.ma[pidx.p, "alt"],
                        gt=rep("HET", length(pidx.p)),
                        sanger=rep("primary", length(pidx.p)),
                        rank=rep(i, length(pidx.p)),
                        is.match.pos=rep(is.match.pos, length(pidx.p)),
                        strand=ps.var.ma[pidx.p, "strand"],
                        trim5=ps.var.ma[pidx.p, "trim5"],
                        trim3=ps.var.ma[pidx.p, "trim3"],
                        query.loc.start=ps.var.ma[pidx.p, "query.loc.start"],
                        query.loc.end=ps.var.ma[pidx.p, "query.loc.end"],
                        query.loc.prime.start=ps.var.ma[pidx.p, "query.loc.prime.start"],
                        query.loc.prime.end=ps.var.ma[pidx.p, "query.loc.prime.end"],
                        query.start=ps.var.ma[pidx.p, "query.start"],
                        query.end=ps.var.ma[pidx.p, "query.end"],
                        ref.prime=ps.var.ma[pidx.p, "ref.prime"],
                        alt.prime=ps.var.ma[pidx.p, "alt.prime"],
                        query.seq.prime=ps.var.ma[pidx.p, "query.seq.prime"],
                        subject.seq.prime=ps.var.ma[pidx.p, "subject.seq.prime"]
                    )                
                }
                if (length(var.pos.ss) > 0) {
                    pidx.s <- which(ss.position %in% var.pos.ss)
                    out.ss <- data.frame(
                        chr=as.character(as.vector(ss.var.ma[pidx.s, "chr"])),
                        start=ss.var.ma[pidx.s, "subject.loc.start"],
                        end=ss.var.ma[pidx.s, "subject.loc.end"],
                        ref=ss.var.ma[pidx.s, "ref"],
                        alt=ss.var.ma[pidx.s, "alt"],
                        gt=rep("HET", length(pidx.s)),
                        sanger=rep("secondary", length(pidx.s)),
                        rank=rep(i, length(pidx.s)),
                        is.match.pos=rep(is.match.pos, length(pidx.s)),
                        strand=ss.var.ma[pidx.s, "strand"],
                        trim5=ss.var.ma[pidx.s, "trim5"],
                        trim3=ss.var.ma[pidx.s, "trim3"],
                        query.loc.start=ss.var.ma[pidx.s, "query.loc.start"],
                        query.loc.end=ss.var.ma[pidx.s, "query.loc.end"],
                        query.loc.prime.start=ss.var.ma[pidx.s, "query.loc.prime.start"],
                        query.loc.prime.end=ss.var.ma[pidx.s, "query.loc.prime.end"],
                        query.start=ss.var.ma[pidx.s, "query.start"],
                        query.end=ss.var.ma[pidx.s, "query.end"],
                        ref.prime=ss.var.ma[pidx.s, "ref.prime"],
                        alt.prime=ss.var.ma[pidx.s, "alt.prime"],
                        query.seq.prime=ss.var.ma[pidx.s, "query.seq.prime"],
                        subject.seq.prime=ss.var.ma[pidx.s, "subject.seq.prime"]
                    )
                }
                out <- list(out.both, out.ps, out.ss)
                out <- out[which(! sapply(out, is.null))]
                if (length(out) == 0) {
                    return(NULL)
                } else {
                    out <- do.call("rbind", out)
                    return(out)
                }
            }
        }
    })
    
    result.lst <- result.lst[which(! sapply(result.lst, is.null))]
    if (length(result.lst) > 0 ) {
        result.lst <- do.call("rbind", result.lst)
    } else {
        result.lst <- NULL
    }
    return(result.lst)
}

## search position
searchPos <- function(blast.dat, pos.lst, sanger.length) {
    check.region <- GRanges(Rle(as.character(pos.lst$chr)), IRanges(start = as.numeric(pos.lst$check.start), end = as.numeric(pos.lst$check.end)), seqinfo = getSeqinfo())
    dat.gr <- blastToGR(blast.dat)
    check.hit <- findOverlaps(dat.gr, check.region)
    if (length(check.hit) == 0) return(NULL)
    idx <- min(queryHits(check.hit))
    ma <- alignMatrix(blast.dat[idx,, drop = FALSE], sanger.length)
    real.pos.gr <- intersect(dat.gr[idx], check.region)
    real.start.index <- which(as.numeric(as.vector(ma[,"subject.loc"])) == start(real.pos.gr))
    real.end.index <- which(as.numeric(as.vector(ma[,"subject.loc"])) == end(real.pos.gr))
    index.vec <- c(real.start.index,real.end.index)
    s.idx <- min(index.vec)
    e.idx <- max(index.vec)
    real.ma <- ma[s.idx:e.idx,, drop=FALSE]
    if (s.idx > e.idx) {
        ref <- paste(rev(real.ma[,"ref"]), collapse = "")
        alt <- paste(rev(real.ma[,"alt"]), collapse = "")
        ref.prime <- paste(rev(real.ma[,"ref.prime"]), collapse = "")
        alt.prime <- paste(rev(real.ma[,"alt.prime"]), collapse = "")
    } else {
        ref <- paste(real.ma[,"ref"], collapse = "")
        alt <- paste(real.ma[,"alt"], collapse = "")
        ref.prime <- paste(real.ma[,"ref.prime"], collapse = "")
        alt.prime <- paste(real.ma[,"alt.prime"], collapse = "")
    }
    out <- data.frame(
                chr=as.character(seqnames(real.pos.gr)),
                start=start(real.pos.gr),
                end=end(real.pos.gr),
                ref=ref,
                alt=alt,
                gt="WT",
                sanger="primary",
                rank=idx,
                is.match.pos=TRUE,
                strand=as.numeric(real.ma[1,"strand"]),
                trim5=min(as.numeric(as.vector(real.ma[,"trim5"]))),
                trim3=min(as.numeric(as.vector(real.ma[,"trim3"]))),
                query.loc.start=min(as.numeric(as.vector(real.ma[,"query.loc"]))),
                query.loc.end=max(as.numeric(as.vector(real.ma[,"query.loc"]))),
                query.loc.prime.start=min(as.numeric(as.vector(real.ma[,"query.loc.prime"]))),
                query.loc.prime.end=max(as.numeric(as.vector(real.ma[,"query.loc.prime"]))),  
                query.start=as.numeric(real.ma[1,"query.start"]),
                query.end=as.numeric(real.ma[1,"query.end"]),
                ref.prime=ref.prime,
                alt.prime=alt.prime,
                query.seq.prime=as.character(real.ma[1,"query.seq.prime"]),
                subject.seq.prime=as.character(real.ma[1,"subject.seq.prime"])   
    )
    return(out)
}

## batch calling variants & search position
batchFunc <- function(hetcall, pos, toolDir, dbPrefix) {
    if (is.null(hetcall)) return(list(het=NULL, variantMatrix=NULL, realp=NULL))
    p.seq <- hetcall@primarySeq
    s.seq <- hetcall@secondarySeq
    seq.len <- length(p.seq)
    p.blast.dat <- doBlast(seq = p.seq, toolDir = toolDir, dbPrefix = dbPrefix)
    s.blast.dat <- doBlast(seq = s.seq, toolDir = toolDir, dbPrefix = dbPrefix)
    var.ma <- NULL
    pos.ma <- NULL
    if.blast <- FALSE
    if (! is.null(p.blast.dat) | ! is.null(s.blast.dat)) {
        if.blast <- TRUE
        if (! is.null(p.blast.dat) & ! is.null(s.blast.dat)) {
            var.ma <- callVariant(ps.dat = p.blast.dat, ss.dat = s.blast.dat, sanger.length = seq.len, pos.lst = pos, trim.bin = 0)
        } else if (is.null(s.blast.dat)) {
            var.ma <- callSingleVariant(var.dat = p.blast.dat, seqtype = "primary", sanger.length = seq.len, pos.lst = pos, trim.bin = 0)
        } else {
            var.ma <- callSingleVariant(var.dat = s.blast.dat, seqtype = "secondary", sanger.length = seq.len, pos.lst = pos, trim.bin = 0)
        }
        if (! is.null(var.ma) & ! is.null(pos)) {
            check.region <- GRanges(Rle(as.character(pos$chr)), IRanges(start = as.numeric(pos$check.start), end = as.numeric(pos$check.end)), seqinfo = getSeqinfo())
            var.ma.gr <- GRanges(Rle(as.character(var.ma$chr)), IRanges(start = as.numeric(var.ma$start), end = as.numeric(var.ma$end)), seqinfo = getSeqinfo())
            hit <- findOverlaps(var.ma.gr, check.region)
            if (length(hit) > 0) pos.ma <- var.ma[queryHits(hit)[1],, drop = FALSE]
        }
        if (is.null(pos.ma) & ! is.null(pos)) pos.ma <- searchPos(blast.dat = p.blast.dat, pos.lst = pos, sanger.length = seq.len)
    }
    return(list(variants=var.ma, positions=pos.ma, if.blast=if.blast))
}

## update hetcalls
updateHetcall <- function(hetcall, winSize, position.df) {
    query.seq <- hetcall@primarySeq
    seq.len <- length(query.seq)
    query.seq.prime.v <- unlist(strsplit(as.character(as.vector(position.df$query.seq.prime)), ""))
    subject.seq.prime.v <- unlist(strsplit(as.character(as.vector(position.df$subject.seq.prime)), ""))    
    query.trim.start <- ifelse (position.df$query.start == 1, "", as.character(query.seq[1:(position.df$query.start-1)]))
    query.trim.end <- ifelse (position.df$query.end == seq.len, "", as.character(query.seq[(position.df$query.end+1):seq.len]))
    query.start.to.loc <- ifelse (position.df$query.loc.prime.start == 1, "", paste(query.seq.prime.v[1:(position.df$query.loc.prime.start-1)], collapse=""))
    query.start.to.loc <- gsub("-", "", query.start.to.loc)
    query.loc.to.end <- ifelse (position.df$query.loc.prime.end == length(query.seq.prime.v), "", paste(query.seq.prime.v[(position.df$query.loc.prime.end+1):length(query.seq.prime.v)], collapse=""))
    query.loc.to.end <- gsub("-", "", query.loc.to.end)
    subject.loc.to.end <- ifelse (position.df$query.loc.prime.start > length(subject.seq.prime.v), "", paste(subject.seq.prime.v[position.df$query.loc.prime.start:length(query.seq.prime.v)], collapse=""))
    # subject.loc.to.end <- gsub("-", "", subject.loc.to.end)
    
    ## new secondary
    new.secondarySeq <- paste0(query.trim.start, query.start.to.loc, gsub("-", "", as.character(as.vector(position.df$alt.prime))), query.loc.to.end, query.trim.end)    
    
    ## new primary
    wins.prime <- ifelse (position.df$query.loc.prime.start - winSize < 1, 1, position.df$query.loc.prime.start - winSize)
    wine.prime <- ifelse (position.df$query.loc.prime.end + winSize > length(query.seq.prime.v), length(query.seq.prime.v), position.df$query.loc.prime.end + winSize)
    query.start.to.win <- ifelse(wins.prime == 1, "", paste(query.seq.prime.v[1:(wins.prime-1)], collapse=""))
    query.win.to.end <- ifelse(wine.prime == length(query.seq.prime.v), "", paste(query.seq.prime.v[(wine.prime+1):length(query.seq.prime.v)], collapse=""))
    query.win.start.to.loc <- ifelse( wins.prime == position.df$query.loc.prime.start, "", paste(query.seq.prime.v[wins.prime:(position.df$query.loc.prime.start-1)], collapse=""))    
    query.win.start.to.loc.clean <- gsub("-", "", query.win.start.to.loc)
    query.loc.to.win.end <- ifelse(wine.prime < position.df$query.loc.prime.start, "", paste(query.seq.prime.v[position.df$query.loc.prime.start:wine.prime], collapse=""))    
    query.loc.to.win.end.clean <- gsub("-", "", query.loc.to.win.end)
    subject.win.start.to.loc <- ifelse( wins.prime == position.df$query.loc.prime.start, "", paste(subject.seq.prime.v[wins.prime:(position.df$query.loc.prime.start-1)], collapse=""))    
    subject.win.start.to.loc.clean <- gsub("-", "", subject.win.start.to.loc)
    subject.loc.to.win.end <- ifelse(wine.prime < position.df$query.loc.prime.start, "", paste(subject.seq.prime.v[position.df$query.loc.prime.start:wine.prime], collapse=""))
    subject.loc.to.win.end.clean <- gsub("-", "", subject.loc.to.win.end)
    
    if (nchar(query.win.start.to.loc) == nchar(query.win.start.to.loc.clean)) {
        cut.1 <- subject.win.start.to.loc
    } else {
        cut.1 <- query.win.start.to.loc.clean
    }
    cut.2 <- ifelse(subject.loc.to.end == "", "", as.character(DNAString(subject.loc.to.end)[1:nchar(query.loc.to.win.end.clean)]))
    new.primarySeq <- paste0(query.trim.start,  gsub("-", "", query.start.to.win), cut.1, cut.2, gsub("-", "", query.win.to.end), query.trim.end)

    new.hetcall <- hetcall
    new.hetcall@primarySeq <- DNAString(new.primarySeq)
    new.hetcall@secondarySeq <- DNAString(new.secondarySeq)
    return(new.hetcall)
}

## functions to plot
seqToLen <- function(file, pos) {
  txt <- pdf_text(file)
  rows <- unlist(strsplit(txt[1], "\n"))
  s <- gsub("^ +", "", rows[1], perl = TRUE)
  ss <- unlist(strsplit(s, ""))
  len <- nchar(s)
  s_count <- 0
  none_count <- 0
  real.pos <- 0
  for (i in 1:length(ss)) {
    if (ss[i] == " ") {
      none_count <- none_count + 1
    } else {
      s_count <- s_count + 1
      if (s_count == pos) real.pos <- i - 1
    }
  }
  return(c(len, real.pos))
}

## add rect
addRect <- function(file, outFile, total.width, bin.width) {  
  img <- image_read_pdf(file, pages = 1, density = 300)
  img_trim <- image_trim(img)
  img_trim.info <- image_info(img_trim)
  img_trim.width <- img_trim.info$width
  img_trim.height <- img_trim.info$height
  crop1 <- image_crop(img_trim, paste0((img_trim.width - 150),"x50+150") )
  crop2 <- image_crop(img_trim, paste0((img_trim.width - 150),"x", (img_trim.height - 50), "+150+50"))
  
  # crop1 & crop2 share the same width
  new_img <- image_append(c(image_border(crop1, "white", "x50"), crop2), stack = TRUE) %>%
    image_annotate("NCBI GRCh37 reference:", location = "+0+10", size = 25, color = "black") %>% 
    image_annotate("Sanger sequence:", location = "+0+120", size = 25, color = "black")
  
  if (bin.width == 10) {
    arrow.pos <- 450
  } else {
    new_infos <- image_info(new_img)
    len_info <- seqToLen(file, pos = bin.width + 1)
    ratio <- len_info[1]/total.width
    arrow.pos <- round(len_info[2] * new_infos$width / (len_info[1] - 1), 0)
  }
  fig <- image_draw(new_img)
  arrows(arrow.pos, 95, arrow.pos, 150, col = "red", lwd = 3)
  dev.off()
  image_write(image_scale(fig, "x750"), outFile, density = 300)
}

## screenshot
screenshotSanger <- function(hetcall, winSize, position.df, outPDF, outPNG) {
    new.hetcall <- updateHetcall(hetcall, winSize = winSize,  position.df = position.df)
    seq.len <- length(new.hetcall@primarySeq)
    realwin <- position.df$query.loc.end - position.df$query.loc.start + 1 + winSize * 2
    plotwin <- 50
    if (winSize > 50) {
        plotwin <- realwin + 1    
    }
    format.trim5 <- position.df$trim5 - winSize
    format.trim3 <- position.df$trim3 - winSize
    if (format.trim5 < 0) format.trim5 <- 0
    if (format.trim3 < 0) format.trim3 <- 0
    chromatogram(new.hetcall, width = plotwin, height = 2, trim5 = format.trim5, trim3 = format.trim3, showcalls = "both", cex.base = 2, cex.mtext = 1, showhets = TRUE, filename = outPDF)
    addRect(outPDF, outPNG, total.width = seq.len, bin.width = winSize)
}

## vcf
makeVcfHeader <- function() {
    fileformat <- "##fileformat=VCFv4.2"
    fileDate <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
    # [END] [GT] [SANGER] [RANK] [IS_MATCH_POS] [STRAND] [TRIM5] [TRIM3] [QUERY_LOC_START] [QUERY_LOC_END] [QUERY_LOC_PRIME_START] [QUERY_LOC_PRIME_END] [QUERY_START] [QUERY_END] [REF_PRIME] [ALT_PRIME] [QUERY_SEQ_PRIME] [SUBJECT_SEQ_PRIME]
    info <- c('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record.">',
    '##INFO=<ID=GT,Number=1,Type=String,Description="HET is short for heterozygous, HOM is short for homozygotes.">',
    '##INFO=<ID=SANGER,Number=1,Type=String,Description="Variant called from sanger primary or secondary sequence, or both sequences.">',
    '##INFO=<ID=RANK,Number=1,Type=Integer,Description="The rank of alignment result to from blast which is used to call this variant, the smaller the better.">',
    '##INFO=<ID=IS_MATCH_POS,Number=1,Type=String,Description="FALSE or TRUE represent if this variant is called from validated region, full stop means unknow.">',
    '##INFO=<ID=STRAND,Number=1,Type=Integer,Description="Strand of sanger sequence.">',
    '##INFO=<ID=TRIM5,Number=1,Type=Integer,Description="The numbers of bases from sanger sequence 5-end to the variant.">',
    '##INFO=<ID=TRIM3,Number=1,Type=Integer,Description="The numbers of bases from the variant to the sanger sequence 3-end.">',
    '##INFO=<ID=QUERY_LOC_START,Number=1,Type=Integer,Description="The start position of variant in query.">',
    '##INFO=<ID=QUERY_LOC_END,Number=1,Type=Integer,Description="The end position of variant in query.">',
    '##INFO=<ID=QUERY_LOC_PRIME_START,Number=1,Type=Integer,Description="The start position of variant in query alignment.">',
    '##INFO=<ID=QUERY_LOC_PRIME_END,Number=1,Type=Integer,Description="The end position of variant in query alignment.">',
    '##INFO=<ID=QUERY_START,Number=1,Type=Integer,Description="The start position of aligned sanger sequence.">',
    '##INFO=<ID=QUERY_END,Number=1,Type=Integer,Description="The end position of aligned sanger sequence.">',
    '##INFO=<ID=REF_PRIME,Number=1,Type=String,Description="The origin variant calling result of reference.">',
    '##INFO=<ID=ALT_PRIME,Number=1,Type=String,Description="The origin variant calling result of sanger sequence.">',
    '##INFO=<ID=QUERY_SEQ_PRIME,Number=1,Type=String,Description="Query alignment blast result.">',
    '##INFO=<ID=SUBJECT_SEQ_PRIME,Number=1,Type=String,Description="Subject alignment blast result.The subject is reference.">')
    
    # GT:GQ:ST:MA:RK:TF:TR
    format <- c('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype.">',
                '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality.">',
                '##FORMAT=<ID=ST,Number=1,Type=Integer,Description="Strand of sanger seq.">',
                '##FORMAT=<ID=MA,Number=1,Type=Integer,Description="0 or 1 represent if this variant is called from validated region, -1 is unknow.">',
                '##FORMAT=<ID=RK,Number=1,Type=Integer,Description="The rank of alignment result to from blast which is used to call this variant, the smaller the better.">',
                '##FORMAT=<ID=TF,Number=1,Type=Integer,Description="The numbers of bases from sanger sequence 5-end to the variant.">',
                '##FORMAT=<ID=TR,Number=1,Type=Integer,Description="The numbers of bases from the variant to the sanger sequence 3-end.">')
    contig <- paste0("##contig=<ID=", seqnames(getSeqinfo()), ",length=", seqlengths(getSeqinfo()), ",assembly=b37>")
    c(fileformat, fileDate, format, info, contig)
}

makeVcfInfo <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { gsub(" ", "", as.character(as.vector(y))) })
    x <- matrix(data = x, ncol = n_col)
    sapply(1:nrow(x), function(i) {
        paste(paste0(toupper(col_n), "=", as.character(as.vector(x[i,]))), collapse=";")
    })
}

makeVcfFormat <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { gsub(" ", "", as.character(as.vector(y))) })
    x <- matrix(data = x, ncol = n_col)
    sapply(1:nrow(x), function(i) {
        paste(as.character(as.vector(x[i,])), collapse=":")
    })
}

# format indel's ref alt 
formatIndel <- function(x) {
    s <- as.numeric(x[1, "start"])
    strand <- as.numeric(x[1, "strand"])
    s.prime <- as.numeric(x[1, "query.loc.prime.start"])
    e.prime <- as.numeric(x[1, "query.loc.prime.end"])
    sub.vec <- DNAString(as.character(x[1, "subject.seq.prime"]))
    ref.minus1 <- as.character(sub.vec[s.prime-1])
    if (strand < 0) {
        ref.minus1 <- as.character(reverseComplement(sub.vec[e.prime+1]))
    }
    alt.list <- strsplit(as.character(x[1, "alt"]), ",")
    ref.vcf <- paste0(ref.minus1, as.character(x[1, "ref"]))
    ref.vcf <- gsub("-", "", ref.vcf)
    alt.vcf <- sapply(alt.list, function(y) { gsub("-", "", paste0(ref.minus1, y)) })
    alt.vcf <- paste(alt.vcf, collapse = ",")
    x[, "start"] <- s - 1
    x[, "ref"] <- ref.vcf
    x[, "alt"] <- alt.vcf
    return(x)
}

makeVcfMatrix <- function(x) {
    CHROM <- as.character(as.vector(x[, "chr"]))
    POS <- as.numeric(as.vector(x[, "start"]))
    ID <- rep(".", nrow(x))
    REF <- gsub("-", ".", as.character(as.vector(x[, "ref"])))
    ALT <- gsub("-", ".", as.character(as.vector(x[, "alt"])))
    REF <- sapply(REF, function(s){
        s <- gsub("B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z", "N", s)
        if (grepl(",", s)) {
            s <- unique(unlist(strsplit(s, ",")))
            return(paste(s, collapse = ","))
        } else {
            return(s)
        }        
    })
    ALT <- sapply(ALT, function(s){
        s <- gsub("B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z", "N", s)
        if (grepl(",", s)) {
            s <- unique(unlist(strsplit(s, ",")))
            return(paste(s, collapse = ","))
        } else {
            return(s)
        }     
    })
    QUAL <- rep(".", nrow(x))
    FILTER <- rep(".", nrow(x))
    info.vec <- setdiff(colnames(x), c("chr", "start", "ref", "alt"))
    new.name <- gsub("\\.", "_", toupper(info.vec))
    info.df <- x[, info.vec, drop = FALSE]
    colnames(info.df) <- new.name
    INFO <- makeVcfInfo(info.df)
    x$gt_format <- ifelse(as.character(as.vector(x[,'gt'])) == "HOM", "1/1", "0/1")
    match_format <- as.character(as.vector(x[,'is.match.pos']))
    match_format[which(match_format == ".")] <- -1
    match_format[which(match_format == "TRUE")] <- 1
    match_format[which(match_format == "FALSE")] <- 0
    x$match_format <- match_format
    x$gq <- rep('.', nrow(x))
    FORMAT <- makeVcfFormat(x[, c("gt_format", "gq", "strand", "match_format", "rank", "trim5", "trim3"), drop = FALSE])
    cbind(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, rep("GT:GQ:ST:MA:RK:TF:TR", nrow(x)), FORMAT)
}

outputVcf <- function(x, id, path) {
    format.x <- lapply(1:nrow(x), function(i) {
        ref <- as.character(x[i,"ref"])
        alt <- as.character(x[i,"alt"])
        alt.list <- unlist(strsplit(alt, ","))
        need.format <- FALSE 
        for (a in alt.list) {
            if (grepl("-", ref) | grepl("-", a) | (nchar(ref) != nchar(a))) need.format <- TRUE
        }
        if (need.format == FALSE) return(x[i,, drop = FALSE])
        return(formatIndel(x[i,, drop = FALSE]))
    })
    format.x <- do.call("rbind", format.x)
    header <- makeVcfHeader()
    vcf.matrix <- makeVcfMatrix(format.x)
    cat(header, file=path, sep="\n")
    cat(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", id), collapse = "\t"),
        file = path, sep = "\n", append = TRUE)
    cat(apply(vcf.matrix, 1, paste, collapse = "\t"), file = path, sep = "\n", append = TRUE)
    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}
