# Make all python utility scripts into executables in the bin directory
# Make all single cell python scripts into executables in the bin directory

# Defile directories
UTILS_DIR := utils
SC_DIR := single_cell_multiome
BIN_DIR := bin
REQUIREMENTS := requirements.txt

# List all python utility scripts
PY_UTILS := $(wildcard $(UTILS_DIR)/*.py)
PY_SC := $(wildcard $(SC_DIR)/*.py)

# List all target executables in the bin directory
EX_UTILS := $(patsubst $(UTILS_DIR)/%.py, $(BIN_DIR)/%, $(PY_UTILS))
EX_SC := $(patsubst $(SC_DIR)/%.py, $(BIN_DIR)/%, $(PY_SC))
EX_ALL := $(EX_UTILS) $(EX_SC)

# Define the default target
all: install $(EX_ALL)

# Define the rule to install dependencies
install: $(REQUIREMENTS)
	pip install -r $(REQUIREMENTS)

# Define the rule to make an executable from a python utility script
$(BIN_DIR)/%: $(UTILS_DIR)/%.py | $(BIN_DIR)
	@echo "Making $@ from $<"

	@# Add shebang, copy to bin directory, and make executable
	@echo "#!/usr/bin/env python" | cat - $< > $@
	@chmod +x $@

# Define the rule to make an executable from a single cell script
$(BIN_DIR)/%: $(SC_DIR)/%.py | $(BIN_DIR)
	@echo "Making $@ from $<"

	@# Add shebang, copy to bin directory, and make executable
	@echo "#!/usr/bin/env python" | cat - $< > $@
	@chmod +x $@

# Create the bin directory if it does not exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Define the clean target
clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean