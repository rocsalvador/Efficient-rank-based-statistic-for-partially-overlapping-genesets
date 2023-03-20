SRC_DIR := src
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cc)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SRC_FILES))

# LDFLAGS := -lpthread
# CCFLAGS := ...
CXXFLAGS := -O3 -Wall
TARGET := gsea
CC := g++

all: obj_dir $(OBJ_FILES)

obj_dir:
	@if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	g++ $(CXXFLAGS) -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG -I/usr/share/R/include -I/home/roc/R/x86_64-pc-linux-gnu-library/4.2/Rcpp/include -fpic  -g -O2 -ffile-prefix-map=/build/r-base-6fXXT3/r-base-4.2.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c -o $@ $<
	g++ -std=gnu++14 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -flto=auto -Wl,-z,relro -o $(TARGET).so $@ -Wl,--export-dynamic -fopenmp -Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -flto=auto -Wl,-z,relro -L/usr/lib/R/lib -lR -lpcre2-8 -llzma -lbz2 -lz -ltirpc -lrt -ldl -lm -licuuc -licui18n -L/usr/lib/R/lib -lR

run:
	Rscript $(SRC_DIR)/$(TARGET).r

clean:
	rm -rf $(OBJ_DIR) $(TARGET).so