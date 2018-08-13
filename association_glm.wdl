task assocTest {

	File gds_file
	File phenotype_file
	String outcome_name
	String covariates_string
	String id_col
	String label
	String test
	File? sample_file
	Int? mac
	String? variant_range

	Int memory
	Int disk

	command {
		echo "Input files" > assocTest.log
		echo "gds_file: ${gds_file}" >> assocTest.log
		echo "phenotype_file: ${phenotype_file}" >> assocTest.log
		echo "outcome_name: ${outcome_name}" >> assocTest.log
		echo "covariates_string: ${covariates_string}" >> assocTest.log
		echo "id_col: ${id_col}" >> assocTest.log
		echo "label: ${label}" >> assocTest.log
		echo "test: ${test}" >> assocTest.log
		echo "sample_file: ${sample_file}" >> assocTest.log
		echo "mac: ${mac}" >> assocTest.log
		echo "variant_range: ${variant_range}" >> assocTest.log
		echo "memory: ${memory}" >> assocTest.log
		echo "disk: ${disk}" >> assocTest.log
		echo "" >> assocTest.log
		dstat -c -d -m --nocolor 10 1>>assocTest.log &
		R --vanilla --args ${gds_file} ${phenotype_file} ${outcome_name} ${covariates_string} ${id_col} ${label} ${test} ${default="NA" sample_file} ${default="5" mac} ${default="NA" variant_range} < /glmAssociation/association_glm.R
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "manninglab/glmassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File assoc = "${label}.assoc.RData"
		File log = "assocTest.log"
	}
}

task summary {
	Float? pval_threshold
	String label
	Array[File] assoc

	Int memory
	Int disk

	command {
		echo "Input files" > summary.log
		echo "pval_threshold: ${pval_threshold}" >> summary.log
		echo "label: ${label}" >> summary.log
		echo "assoc: ${sep = ',' assoc}" >> summary.log
		echo "memory: ${memory}" >> summary.log
		echo "disk: ${disk}" >> summary.log
		echo "" >> summary.log
		dstat -c -d -m --nocolor 10 1>>summary.log &
		R --vanilla --args ${default="0.0001" pval_threshold} ${label} ${sep = ',' assoc} < /glmAssociation/summary.R
	}
	
	runtime {
		docker: "manninglab/glmassociation:latest"
  	    disks: "local-disk ${disk} SSD"
        memory: "${memory} GB"
	}

	output {
		File allassoccsv = "${label}.assoc.csv"
		File topassoccsv = "${label}.topassoc.csv"
		File plots = "${label}_association_plots.png"
		File log = "summary.log"
	}
}

workflow w_assocTest {
	Array[File] these_gds_file
	File this_phenotype_file
	String this_outcome_name
	String this_covariates_string
	String this_id_col
	String this_label
	String this_test
	File? this_sample_file
	Int? this_mac
	String? this_variant_range

	Float? this_pval_threshold

	Int this_memory
	Int this_disk

	scatter(this_gds_file in these_gds_file) {


		call assocTest {
			input: gds_file = this_gds_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, covariates_string = this_covariates_string, id_col = this_id_col, label = this_label, test = this_test, sample_file = this_sample_file, mac = this_mac, variant_range = this_variant_range, memory = this_memory, disk = this_disk
		}

	}
	
	call summary {
		input: pval_threshold = this_pval_threshold, label = this_label, assoc = assocTest.assoc, memory = this_memory, disk = this_disk
	}
}