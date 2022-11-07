#!/bin/bash

## F19K16_F24B22only
#bash ./get_mock_fastqs_wBAC.sh \
#-d "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data" \
#-f "/cluster/home/t110409uhn/git/make_toy_fastqs/backups/fastq_flist.bak" \
#-b "/cluster/home/t110409uhn/git/make_toy_fastqs/backups/bedpe_flist.bak" \
#-m "toy_sample1_F19K16_F24B22only" \
#-i "ABCD:ABCD" \
#-n 1000000 \
#-s 42 \
#-c "BACs" \
#-o "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_BACs/toy_sample1" \
#-t
## { ~400k reads from F19K16_F24B22 }

## hg38_F19K16_F24B22
bash ./get_mock_fastqs_wBAC.sh \
-d "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data" \
-f "/cluster/home/t110409uhn/git/make_toy_fastqs/backups/fastq_flist.bak" \
-b "/cluster/home/t110409uhn/git/make_toy_fastqs/backups/bedpe_flist.bak" \
-m "toy_sample1" \
-i "AAAA:AAAA" \
-n 400000 \
-s 111 \
-c "all" \
-o "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_allchrs/toy_sample1" \
-t

#bash ./get_mock_fastqs_wBAC.sh \
#-d "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data" \
#-f "/cluster/home/t110409uhn/git/make_toy_fastqs/backups/fastq_flist.bak" \
#-b "/cluster/home/t110409uhn/git/make_toy_fastqs/backups/bedpe_flist.bak" \
#-m "toy_sample2" \
#-i "BBBB:BBBB" \
#-n 400000 \
#-s 222 \
#-c "all" \
#-o "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_allchrs/toy_sample2" \
#-t
