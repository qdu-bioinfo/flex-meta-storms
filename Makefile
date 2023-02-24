CC:=g++

ifneq ($(shell uname -s),Darwin)
    # Not macOS
else ifneq ($(shell which g++-12),)
    # g++-12 exists
    CC:=g++-12
else ifneq ($(shell which g++-11),)
    # g++-11 exists
    CC:=g++-11
else ifneq ($(shell which g++-10),)
    # g++-10 exists
    CC:=g++-10
else ifneq ($(shell which g++-9),)
    # g++-9 exists
    CC:=g++-9
else
    # g++-12, g++-11, g++-10, and g++-9 do not exist
endif

OMPFLG=-fopenmp
HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched
EXE_CMP=bin/FMS-comp-taxa

tax:$(OBJ_TAX) src/fms_comp_sam.cpp
	$(CC) -o $(EXE_CMP) src/fms_comp_sam.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	chmod +x Rscript/PM_Marker_Test.R
clean:
	rm -rf bin/PM-* src/*.o
