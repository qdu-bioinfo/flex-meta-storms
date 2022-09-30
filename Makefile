CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
	exist = $(shell if [ -e '/usr/local/bin/g++-10' ]; then echo "exist"; else echo "notexist"; fi;)
	ifeq ($(exist),exist)
		CC:=g++-10
	else
        	exist = $(shell if [ -e '/usr/local/bin/g++-9' ]; then echo "exist"; else echo "notexist"; fi;)
        	ifeq ($(exist),exist)
                	CC:=g++-9
		else
			CC:=g++-8
		endif
	endif
endif
OMPFLG=-fopenmp
HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse
EXE_CMP=bin/FMS-comp-taxa

tax:$(OBJ_TAX) src/fms_comp_sam.cpp src/comp_sam.cpp
	$(CC) -o $(EXE_CMP) src/fms_comp_sam.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
clean:
	rm -rf bin/PM-* src/*.o
