CCDIR := src/cc
RCPPDIR := src/rcpp

all:
	@make -C $(CCDIR)
	@make -C $(RCPPDIR)

cc:
	@make -C $(CCDIR)

runcc:
	@$(CCDIR)/gsea

rcpp:
	@make -C $(RCPPDIR)

clean:
	@make clean -C $(CCDIR)
	@make clean -C $(RCPPDIR)