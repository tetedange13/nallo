include { STRAGLR          } from '../../../modules/local/straglr/'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main/'
include { ADD_FOUND_IN_TAG } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_MERGE   } from '../../../modules/nf-core/bcftools/merge/'

workflow CALL_REPEAT_EXPANSIONS_STRAGLR {

    take:
    ch_bam_bai  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_bed      // channel: [mandatory] [ val(meta), path(bed) ]

    main:
    ch_versions = Channel.empty()

    STRAGLR (
        ch_bam_bai,
        ch_fasta,
        ch_bed.map { _meta, bed -> bed }
    )
    ch_versions.mix(STRAGLR.out.versions)

    TABIX_BGZIPTABIX(STRAGLR.out.vcf)
    ch_vcf = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, _tbi -> [ meta, vcf ] }
    ch_tbi = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _vcf, tbi -> [ meta, tbi ] }

    ADD_FOUND_IN_TAG (
        TABIX_BGZIPTABIX.out.gz_tbi,
        "STRaglr"
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    ADD_FOUND_IN_TAG.out.vcf
        .join(ADD_FOUND_IN_TAG.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [ [ id: meta.family_id ], vcf, tbi ] }
        .groupTuple()
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE (
        ch_bcftools_merge_in,
        [ [], [] ],
        [ [], [] ],
        [ [], [] ]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    sample_vcf  = ch_vcf                   // channel: [ val(meta), path(vcf) ]
    sample_tbi  = ch_tbi                   // channel: [ val(meta), path(tbi) ]
    family_vcf  = BCFTOOLS_MERGE.out.vcf   // channel: [ val(meta), path(vcf) ]
    family_tbi  = BCFTOOLS_MERGE.out.index // channel: [ val(meta), path(tbi) ]
    versions    = ch_versions              // channel: [ versions.yml ]

}
