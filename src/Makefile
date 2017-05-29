# Copyright (C) 2015 Andrew D. Smith
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

all_subdirs=prog common
app_subdirs=prog 

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(realpath smithlab_cpp)
endif

ifndef TREETOOL
TREETOOL=$(realpath adssrc/treetool)
endif

ifndef BOOST
BOOST=$(realpath boost_1_64_0)
endif

all:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) BOOST=$(BOOST) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) BOOST=$(BOOST) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=${TREETOOL} BOOST=$(BOOST) test; \
	done;

clean:
	@for i in $(all_subdirs); do \
	make -C $${i} clean; \
	done;

.PHONY: test

