#### gtf processing 
library(Rsubread)
library(rtracklayer)


mice_gtf <- rtracklayer::import('Mus_musculus.GRCm39.110.gtf')

mice_gtf <- rtracklayer::as.data.frame(mice_gtf)

mice_gtf_canonicals <- mice_gtf[which(mice_gtf$tag == "Ensembl_canonical"),]


mice_gtf_canonicals_filtered_500 <- Reduce(rbind,
                                           lapply(unique(mice_gtf_canonicals$transcript_id), function(x) {
                                             transcript_ranges <- mice_gtf_canonicals[which(mice_gtf_canonicals$transcript_id == x),]
                                             
                                             transcript_exones <- transcript_ranges[which(transcript_ranges$type == "exon"),]
                                             transcript_exones <- transcript_exones[order(as.integer(transcript_exones$exon_number), decreasing = T),]
                                             transcript <- transcript_ranges[which(transcript_ranges$type == "transcript"),]
                                             
                                             rem_len <- 500
                                             if(transcript$strand == "+") {
                                               Reduce(rbind,
                                                      lapply(1:nrow(transcript_exones), function (y) {
                                                        if(rem_len > 0) {
                                                          transcript_exones[y,"start"] <- transcript_exones[y,"end"] - min(transcript_exones[y,"width"] - 1, rem_len)
                                                          rem_len <<- rem_len - min(transcript_exones[y,"width"] - 1, rem_len)
                                                          transcript_exones[y,]
                                                        }
                                                      })
                                               )
                                             } else {
                                               Reduce(rbind,
                                                      lapply(1:nrow(transcript_exones), function (y) {
                                                        if(rem_len > 0) {
                                                          transcript_exones[y,"end"] <- transcript_exones[y,"start"] + min(transcript_exones[y,"width"] - 1, rem_len)
                                                          rem_len <<- rem_len - min(transcript_exones[y,"width"] - 1, rem_len)
                                                          transcript_exones[y,]
                                                        }
                                                      })
                                               )
                                             }
                                           })
)

rtracklayer::export(mice_gtf_canonicals_filtered_500, "3_prime_filtered_gtf_canonical_transcripts_500_fixed_exon_order.gtf", format = "gtf")