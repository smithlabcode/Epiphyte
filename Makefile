# Copyright (C) 2015 University of Southern California and
#                    Andrew D. Smith
#
# Authors: Jenny Qu, Andrew D. Smith
#
# This code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#


EPIPHYTE_ROOT = $(shell pwd)
export PATH := $(shell pwd):$(PATH)
BINDIR = $(EPIPHYTE_ROOT)/bin

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/src/smithlab_cpp
ifeq ("$(wildcard $(SMITHLAB_CPP))","")
$(error SMITHLAB_CPP variable not set and smithlab_cpp not found)
endif
endif

ifndef TREETOOL
TREETOOL=$(abspath $(dir $(MAKEFILE_LIST)))/src/adssrc/treetool
ifeq ("$(wildcard $(TREETOOL))","")
$(error could not set TREETOOL variable)
endif
endif

all:
	@make -C src EPIPHYTE_ROOT=$(EPIPHYTE_ROOT) OPT=1

install:
	@export PATH=$(PATH)
	@make -C src EPIPHYTE_ROOT=$(EPIPHYTE_ROOT) OPT=1 install

test:
	@make -C src OPT=1 test
.PHONY: test

clean:
	@rm -rf $(BINDIR)
	@make -C src EPIPHYTE_ROOT=$(EPIPHYTE_ROOT) clean
.PHONY: clean
