SRC_DIR := .
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cc)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SRC_FILES))

LDFLAGS := -lpthread
# CCFLAGS := ...
CXXFLAGS := -O3 -Wall
TARGET := gsea
CC := g++

all: obj_dir $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CC) $(LDFLAGS) -o $@ $^

obj_dir:
	@if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CC) $(CCFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
