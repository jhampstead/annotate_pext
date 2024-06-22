# annotate_pext
Fast proportion expressed across transcripts (PEXT) score calculation from VEP-annotated VCF files.

The calculation of PEXT scores is based on the 2020 paper *[Transcript expression-aware annotation improves rare variant interpretation](https://www.nature.com/articles/s41586-020-2329-2)* by Cummings *et al.*

## Running annotate_pext

```annotate_pext``` is written in C, and uses HTSLib functions. HTSLib is required to compile and execute.

If an HTSLib module is already available on your system, it should be loaded before compiling or running the executable:

```
module spider htslib
module load bioinf/htslib/1.9
```

If HTSLib is unavailable on your system, you can [download HTSLib here](https://www.htslib.org/download/) and [install from source](https://github.com/samtools/htslib/blob/develop/INSTALL). If ```extract_vcf``` can't find HTSLib at runtime, execution will fail with ```./annotate_pext: error while loading shared libraries: libhts.so.2: cannot open shared object file: No such file or directory```. I compile and run with HTSLib 1.9; other versions of HTSLib have not been tested and may cause unexpected behaviour.

```annotate_pext``` requires the following command line arguments:
* An isoform expression matrix <isoform_expression_matrix.tsv>. The rows of this file should represent transcripts, while the columns should represent tissues. This script has been optimised for use with the GTEx isoform expression matrix provided at ```test/GTEx_median_tissue_expression_matrix.tsv```. Other expression matrices should adhere to this format.
* An input VCF, or one piped through STDIN
* An output VCF, or redirection through STDOUT

Both input and output VCFs can be uncompressed or uncompressed, indexed or unindexed.

By default PEXT scores are calculated for all coding variants (defined as variants within the CDS of each transcript) across all tissues in the provided isoform expression matrix. The following optional arguments modify this behaviour:
* ```--most severe <ordered_consequences.txt>``` calculates PEXT scores only for the most severe consequence per variant. For example, if a variant has synonymous_variant, synonymous_variant, synonymous_variant, stop_gained, stop_gained consequences across the five transcripts covering the position, PEXT scores would only be calculated for stop_gained even though synonymous_variants are also coding.
* ```--tissue_group [TISSUE_GROUP]``` calculates PEXT scores only across certain tissues within the isoform expression matrix. To calculate PEXT across only Brain tissues, for example, ```--tissue-group Brain``` would use only Brain - Amygdala, Brain - Anterior cingulate cortex (BA24), Brain - Caudate (basal ganglia), Brain - Cerebellar Hemisphere, Brain - Cerebellum, Brain - Cortex, Brain - Frontal Cortex (BA9), Brain - Hippocampus, Brain - Hypothalamus, Brain - Nucleus accumbens (basal ganglia), Brain - Putamen (basal ganglia), Brain - Spinal cord (cervical c-1), and Brain - Substantia nigra tissues.

To test whether ```annotate_pext``` is working correctly, you can run the following command:

```
./annotate_pext test/GTEx_median_tissue_expression_matrix.tsv test/test.vcf test/test_out.vcf
```

## Compiling annotate_pext

If you need to compile ```annotate_pext``` from source, first load HTSLib. You can then compile the executable using the following command:

```gcc -O3 -o annotate_pext annotate_pext.c -lhts```

Compilation was done using gcc 10.2.0 on Linux.

