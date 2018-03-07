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
		R --vanilla --args ${gds_file} ${phenotype_file} ${outcome_name} ${covariates_string} ${id_col} ${label} ${test} ${default="NA" sample_file} ${default="5" mac} ${default="NA" variant_range} < /glmAssociation/association_glm.R
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "manning-lab/glmAssociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
	String test
	Float? pval_threshold
	String label
	Array[File] assoc

	Int memory
	Int disk

	command {
		R --vanilla --args ${test} ${default="0.0001" pval_threshold} ${label} ${sep = ',' assoc} < /glmAssociation/summary.R
	}
	
	runtime {
		docker: "manning-lab/glmAssociation:latest"
  	    disks: "local-disk ${disk} SSD"
        memory: "${memory}G"
	}

	output {
		File allassoccsv = "${label}.assoc.csv"
		File topassoccsv = "${label}.topassoc.csv"
		File plots = "${label}_association_plots.png"
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
		input: test = this_test, pval_threshold = this_pval_threshold, label = this_label, assoc = assocTest.assoc, memory = this_memory, disk = this_disk
	}
