DIR_INC = ./include
DIR_SRC = ./source
DIR_OBJ = ./obj


SRC = $(wildcard ${DIR_SRC}/*.cu)\
		$(wildcard ${DIR_SRC}/*.cpp)
NO_DIR = $(notdir $(SRC))
OBJ := $(patsubst %.cu, $(DIR_OBJ)/%.o,$(NO_DIR))
OBJ := $(patsubst %.cpp, $(DIR_OBJ)/%.o, $(OBJ))

TARGET = TwoElec.out


NVCC = nvcc
NVCC_FLAGS =  -Xcompiler -Wall -O4  -lcurand -rdc=true -I $(DIR_INC) 




$(sort $(DIR_OBJ)) :
	mkdir -p $@

#test:
#	 @echo $(SRC)
#	 @echo $(NO_DIR)
#	 @echo $(OBJ)
#	 @echo "end"


${TARGET}:${OBJ}
	$(NVCC) $(OBJ)  -o $@
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.cu
	$(NVCC) $(NVCC_FLAGS) -c  $< -o $@
${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	$(NVCC) $(NVCC_FLAGS) -c  $< -o $@
	
.PHONY:clean
clean:
	@find ${DIR_OBJ} -name *.o -exec rm -rf {} \;
	@find ${TARGET} -exec rm -rf {} \;
