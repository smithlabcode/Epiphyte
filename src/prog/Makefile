#    Copyright (C) 2015 Andrew D. Smith
#
#    Authors: Jenny Qu, Andrew D. Smith
#
#    This code is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(realpath ../smithlab_cpp)
endif

ifndef TREETOOL
TREETOOL=$(realpath ../adssrc/treetool)
endif

ifndef METHPIPE_DIR
METHPIPE_DIR=$(METHPIPE_ROOT)/src/common
endif


ifndef BOOST_ROOT
ifeq ("$(wildcard /usr/include/boost)","")
$(error need to set variable BOOST_ROOT)
else
BOOST_ROOT = /usr/include/
endif	
endif 

PROGS = indep-methprob indep-epi-phylo epiphy-sim epiphy-est \
	epiphy-post epiphy-seg 

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50 -std=c++11
OPTFLAGS = -O2
DEBUGFLAGS = -g
LIBS = -lgsl -lgslcblas -L/usr/local/lib 

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = $(realpath ../common)
INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR) ${TREETOOL} $(METHPIPE_DIR) ${BOOST_ROOT} /usr/local/include
INCLUDEARGS = $(addprefix -I, $(INCLUDEDIRS))

all: $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

indep-epi-phylo: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o \
	BetaBin.o PhyloTreePreorder.o epiphy_utils.o ) \
	$(addprefix $(TREETOOL)/, PhyloTree.o ) \
	$(addprefix $(METHPIPE_DIR)/, MethpipeSite.o)

epiphy-est epiphy-post epiphy-seg epiphy-est-multichain: \
	$(addprefix $(COMMON_DIR)/, MethpipeFiles.o \
	param_set.o PhyloTreePreorder.o sufficient_statistics_helpers.o \
	optimize_params.o epiphy_utils.o epiphy_mcmc.o) \
	$(addprefix $(METHPIPE_DIR)/, MethpipeSite.o) \
        $(addprefix $(TREETOOL)/, PhyloTree.o)

epiphy-sim: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o PhyloTreePreorder.o \
	param_set.o PhyloTreePreorder.o sufficient_statistics_helpers.o \
	optimize_params.o) \
	$(addprefix $(METHPIPE_DIR)/, MethpipeSite.o) \
        $(addprefix $(TREETOOL)/, PhyloTree.o)

indep-methprob: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o BetaBin.o)

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
