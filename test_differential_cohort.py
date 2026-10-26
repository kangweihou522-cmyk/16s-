"""Real-R cohort and descriptive-abundance checks; no ANCOM model is fitted."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
RSCRIPT = os.environ.get('WORKSTATION_TEST_RSCRIPT') or shutil.which('Rscript')
if not RSCRIPT and Path('/workspace/r-conda-review/Rscript').is_file():
    RSCRIPT = '/workspace/r-conda-review/Rscript'


@unittest.skipUnless(RSCRIPT, 'A real R runtime is not installed')
class DifferentialCohortTests(unittest.TestCase):
    def run_r(self, expression):
        with tempfile.TemporaryDirectory(prefix='16s-diff-cohort-') as output:
            result = subprocess.run([RSCRIPT, '-e', expression,
                                     str(ROOT / 'differential_ancombc.R'), output],
                                    capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_library_boundary_prevalence_order_and_insufficient_groups(self):
        self.run_r("""a<-commandArgs(TRUE);source(a[1])
ids<-c('0001','0002','low_A','B1','B2','B3')
m<-rbind(common_1=c(6,8,2,12,20,10),common_2=c(4,2,6,8,0,0),rare=c(1,0,1000,0,0,0));colnames(m)<-ids
meta<-data.frame(Group=factor(rep(c('A','B'),each=3),levels=c('A','B')),row.names=ids)
o<-list(prv_cut=.5,lib_cut=10)
c<-ancom_sample_cohort(m,meta,o);validate_ancom_sample_cohort(c)
stopifnot(identical(colnames(c$counts),ids[c(1,2,4,5,6)]),identical(rownames(c$meta),colnames(c$counts)),
          c$audit$input_library_size[3]==1008,c$audit$library_size_after_prevalence[3]==8,
          c$audit$reason[3]=='below_lib_cut_after_prevalence',
          all(c$audit$retained_by_library_filter[c$audit$library_size_after_prevalence==10]),
          c$summary$status=='planned_before_model',is.null(c$summary$model_sample_ids))
fit<-list(feature_table=m[1:2,colnames(c$counts),drop=FALSE]);actual<-ancom_sample_cohort(m,meta,o,fit)
stopifnot(actual$summary$source=='official_fit_feature_table',actual$summary$model_sample_count==5)
write_ancom_sample_cohort(actual,a[2]);saved<-read.delim(file.path(a[2],'sample_inclusion.tsv'),colClasses='character')
stopifnot(identical(saved$SampleID,ids),saved$reason[3]=='below_lib_cut_after_prevalence')
failure<-tryCatch(validate_ancom_sample_cohort(ancom_sample_cohort(m,meta,modifyList(o,list(lib_cut=11)))),error=conditionMessage)
stopifnot(grepl('at least two',failure),grepl('A=0',failure))
bad<-tryCatch(ancom_sample_cohort(m,meta,o,list(feature_table=m)),error=conditionMessage)
stopifnot(grepl('differ',bad))
missing<-tryCatch(ancom_sample_cohort(m,meta,o,list(res=data.frame())),error=conditionMessage)
stopifnot(grepl('feature_table is required',missing))
default<-ancom_sample_cohort(m,meta,modifyList(o,list(lib_cut=0)),list(res=data.frame()))
stopifnot(identical(default$counts,m),identical(default$meta,meta))
cat('Cutoff boundary, prevalence ordering, exact samples and group guard verified without fitting a model')""")

    def test_abundance_uses_model_samples_with_full_original_taxon_denominators(self):
        self.run_r("""a<-commandArgs(TRUE);source(a[1])
ids<-c('0001','0002','low_A','B1','B2','B3')
m<-rbind(common_1=c(6,8,2,12,20,10),common_2=c(4,2,6,8,0,0),rare=c(1,0,1000,0,0,0));colnames(m)<-ids
meta<-data.frame(Group=factor(rep(c('A','B'),each=3),levels=c('A','B')),row.names=ids)
res<-data.frame(taxon='common_1',lfc_GroupB=1,se_GroupB=.1,p_GroupB=.001,q_GroupB=.01,
                diff_GroupB=TRUE,passed_ss_GroupB=TRUE,check.names=FALSE)
fit<-list(res=res,feature_table=m[1:2,c(1,2,4,5,6),drop=FALSE])
o<-list(prv_cut=.5,lib_cut=10,reference_group='A',alpha=.05,pseudo_sens=TRUE,font='sans')
report<-plot_ancom_level(fit,m,meta,a[2],o,'Synthetic saved model')
means<-read.delim(file.path(a[2],'abundance','group_mean_relative_abundance.tsv'),check.names=FALSE)
stopifnot(abs(means$A[means$taxon=='common_1']-mean(c(6/11,8/10)))<1e-12,
          report$sample_cohort$model_sample_count==5,report$sample_cohort$group_sizes$A==2)
sample<-read.delim(file.path(a[2],'abundance','sample_relative_abundance.tsv'),check.names=FALSE)
stopifnot(!'low_A'%in%names(sample),'rare'%in%sample$taxon,abs(sum(sample[['0001']])-1)<1e-12)
stopifnot(grepl('retained ANCOM-BC2 model sample',paste(readLines(file.path(a[2],'abundance','abundance_heatmap.svg')),collapse=' ')))
cat('Model cohort plotting retains full original taxon denominators; no model fitted')""")


if __name__ == '__main__':
    unittest.main()
