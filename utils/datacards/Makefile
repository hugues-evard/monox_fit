# Determine all needed paths
MAKEFILE_PATH := $(realpath $(lastword $(MAKEFILE_LIST)))
TOPDIR := $(realpath $(dir $(MAKEFILE_PATH))/..)
YEAR := $(notdir $(abspath $(CURDIR)/../../))
CHANNEL := $(notdir $(abspath $(CURDIR)/../../../))
SCRIPTDIR := $(TOPDIR)/datacards

# List of supported script targets
SCRIPTS := diagnostics impacts impacts_condor limits gof

# Declare all phony targets
.PHONY: $(SCRIPTS) plots printcommand-%

# Generic rule to run scripts matching the target name: do_<target>.sh
$(SCRIPTS):
	@echo "Running script for $@ on channel=$(CHANNEL)"
	bash $(SCRIPTDIR)/do_$@.sh $(CHANNEL)

# Rule to print the content of a script
printcommand-%:
	@echo "Contents of do_$*.sh:"
	cat $(SCRIPTDIR)/do_$*.sh
	@echo ""
	@echo "End of do_$*.sh"

# Alias for diagnostics
.PHONY: diag printcommand-diag
diag: diagnostics
printcommand-diag: printcommand-diagnostics

# Create datacards
cards:
	@echo "Generating cards for channel=$(CHANNEL) year=$(YEAR)"
	python3 $(SCRIPTDIR)/build_cards.py --channel $(CHANNEL) --year $(YEAR)

# Plotting
plots:
	@echo "Generating plots for channel=$(CHANNEL) year=$(YEAR)"
	python3 $(SCRIPTDIR)/make_plots.py --channel $(CHANNEL) --year $(YEAR)
