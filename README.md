# mlift 
A mummer based genomic coordinate liftover tool

This is a working but prototype version so use this with caution. If you provide a genome coordinate bed file (or simply a tab delimited file), mlift will lift the coordinates of ref genome to query genome (ref and query genome assignments are arbitrary and follow the MUMmer convention). If you use it, and has a feedback, please email me at mchakrab@uci.edu.

1. Install

 ```
	Make
 ```
2. Run nucmer

The parameters of nucmer depends on how you want to lift your coordinates. If you want to ignore duplicates, use 

  ```
	nucmer -mum -prefix ref2q ref.fasta query.fasta

	delta-filter -r -q re2q.delta >ref2q.rq.delta
 ```
If you want to get all copies, use

 ```
	nucmer -prefix ref2q.mm ref.fasta query.fasta
 ```

3. Run mlift

 ```
	mlift ref2q.delta foo.bed s
 ```
Here the delta file is the delta alignment from the nucmer run in step 2. You can use 's' or 'l' as the mode. When your complete interval is not present in the alignment, 's' will contract the interval to fit it inside the alignment, whreas 'l' will write 'NA' for your coordinate that is not covered in the alignment. I use 's' when I have to find partial of full gene duplications using mlift. The bed file should have the ref genome coordinates in the following tab separated format -
	
	Chromosome_name	Start	End
