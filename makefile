SRC_DIR := src
TARGET := gseacc

build:
	rm -rf $(TARGET)
	R -e 'library("Rcpp");filenames <- c(Sys.glob("$(SRC_DIR)/*.cc"), Sys.glob("$(SRC_DIR)/*.hh"));Rcpp.package.skeleton("$(TARGET)", cpp_files = filenames)'
	R CMD build $(TARGET)

install: build
	R CMD INSTALL $(TARGET)_1.0.tar.gz 

run:
	Rscript $(SRC_DIR)/gsea.r

clean:
	rm -rf $(TARGET)