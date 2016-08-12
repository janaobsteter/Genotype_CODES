#############################################################
#script to extract microsatelites from phased Beagle files
#run from Beagle_imputedGP4 directory
###############################################################

	awk '$2== "BM1824" || $2=="id" {print}' chr1.Chr1.bgl.phased > 1a.txt
	awk '$2== "BM2113" || $2=="id" {print}' chr2.Chr2.bgl.phased > 2a.txt
	awk '$2== "INRA023" || $2=="id" {print}' chr3.Chr3.bgl.phased > 3a.txt
	awk '$2== "ETH10" || $2=="id" {print}' chr5.Chr5.bgl.phased > 5a.txt
	awk '$2== "ETH225" || $2=="id" {print}' chr9.Chr9.bgl.phased > 9a.txt
	awk '$2== "SPS115" || $2=="id" {print}' chr15.Chr15.bgl.phased > 15a.txt
	awk '$2== "TGLA53" || $2=="id" {print}' chr16.Chr16.bgl.phased > 16a.txt
	awk '$2== "TGLA227" || $2=="id" {print}' chr18.Chr18.bgl.phased > 18a.txt
	awk '$2== "ETH3" || $2=="id" {print}' chr19.Chr19.bgl.phased > 19a.txt
	awk '$2== "TGLA126" || $2=="id" {print}' chr20.Chr20.bgl.phased > 20a.txt
	awk '$2== "TGLA122" || $2=="id" {print}' chr21.Chr21.bgl.phased > 21a.txt
	awk '$2== "BM1818" || $2=="id" {print}' chr23.Chr23.bgl.phased > 23a.txt



