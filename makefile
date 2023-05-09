SRC_DIR := src
R_DIR := R
TARGET := gseacc

cc:
	g++ -O3 -Wall -c src/main.cc src/gsea.cc
	g++ -o gsea main.o gsea.o

build:
	rm -rf $(TARGET)
	R -e 'library("Rcpp");filenames <- c(Sys.glob("$(SRC_DIR)/*.cc"), Sys.glob("$(SRC_DIR)/*.hh"));Rcpp.package.skeleton("$(TARGET)", cpp_files = filenames, code_files = "$(R_DIR)/gseacc.R")'
	R CMD build $(TARGET)

install: build
	R CMD INSTALL $(TARGET)_1.0.tar.gz 

run:
	Rscript $(SRC_DIR)/gsea.r &> log-$(shell date +%s).txt

run-sc:
	Rscript src/scrna-gsea.r &> log.txt

debug:
	g++ -g -c src/main.cc src/gsea.cc
	g++ -o gsea main.o gsea.o

clean:
	rm -rf $(TARGET) $(TARGET)_1.0.tar.gz *.o gsea
