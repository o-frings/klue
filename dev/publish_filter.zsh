# dev/publish_filter.zsh -- sourced by dev/make_reproduction_archive.sh: the
# result files the archive publishes, the same set this repository tracks.
#
# keep_output <file name>: 0 = publish this top-level output/ file.
#   not published:
#   *.log *.done                 run logs and completion markers
#   *_chunk[0-9][0-9].rds        resume part-files of the checkpointed runners
#   compare_klue_versions_*      the klue-version equivalence harness (~10 MB)
#   benchmark_availability* benchmark_concomitant* benchmark_lcmnl_se*
#   validate_* diagnose_* compare_packages_stress.rds
#                                development cross-checks against Apollo and
#                                mlogit, not paper results
#   validation_* smoketest_*     adapter self-test artifacts, not paper results
#   standard_rerun_summary*      revision bookkeeping against the pre-revision
#                                version12.tex; it needs that manuscript, the lane
#                                logs and output/superseded_klue095/
#   vittel_respec_lcmnl.rds      the Vittel ladder's fit objects, which carry each
#                                respondent's class posteriors (930 x C per fit);
#                                respondent-level Vittel results are not
#                                redistributed, as test_data_files/ is not. The
#                                table is published as vittel_respec_lcmnl.csv.
#   (output/benchmark_apollo_defaults.* is a paper result and is published, and so
#   are the two stress-ladder logs, the record of what Apollo printed in each cell:
#   the apollo_lcEM iteration cap and the numerical-gradient message.)
# OUTPUT_SUBDIRS: output/ subfolders published whole (their .rds and .csv):
#   mmnl_bench_std/  the empirical MMNL fits behind the benchmark tables
#   mmnl_bench/      the klue 0.9.5 fits that README_REPRODUCE.md still cites
keep_output() {
  case "$1" in
    stress_replicate_blocked_std.log|stress_replicate_std.log) return 0 ;;
    *.log|*.done|v10_backup) return 1 ;;
    *_chunk[0-9][0-9].rds) return 1 ;;
    compare_klue_versions_*|compare_packages_stress.rds) return 1 ;;
    benchmark_availability*|benchmark_concomitant*|benchmark_lcmnl_se*) return 1 ;;
    validate_*|diagnose_*) return 1 ;;
    validation_*|smoketest_*) return 1 ;;
    standard_rerun_summary*) return 1 ;;
    vittel_respec_lcmnl.rds) return 1 ;;
    *.csv|*.txt|*.rds) return 0 ;;
    *) return 1 ;;
  esac
}
OUTPUT_SUBDIRS=(mmnl_bench_std mmnl_bench)
