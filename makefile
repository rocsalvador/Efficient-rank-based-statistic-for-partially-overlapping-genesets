SRC_DIR := src
TARGET := gsea

install:
	@rm -rf gseacc
	@R -e 'library("Rcpp");filenames <- c(Sys.glob("src/*.cc"), Sys.glob("src/*.hh"));Rcpp.package.skeleton("gseacc", cpp_files = filenames)'
	@R CMD build gseacc
	@R CMD INSTALL gseacc_1.0.tar.gz 

run:
	Rscript $(SRC_DIR)/$(TARGET).r

clean:
	rm -rf $(OBJ_DIR) $(TARGET).so